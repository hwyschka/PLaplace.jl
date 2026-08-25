"""
    mutable struct IterationData
    
Stores all data that changes during the iteration of the interior-point method.
Especially computes gradient and hessian of the barrier function in each iteration.
"""
mutable struct IterationData
    "iteration vector"
    x::AbstractVector{Float64}
    "current descent direction"
    dx::AbstractVector{Float64}

    "barrier function"
    barrier::BarrierFunction
    "gradient of the barrier function at x"
    gradientF::AbstractVector{Float64}
    "factorized hessian of the barrier function at x"
    hessianF::Any
    "preconditioner for the hessian"
    P::Any
end

function IterationData(x::AbstractVector{Float64}, bf::BarrierFunction)
    le = length(x)
    IterationData(
        x,
        zeros(Float64, le),
        bf,
        zeros(Float64, le),
        spzeros(Float64, le, le),
        LinearAlgebra.I
    )
end

function IterationData(bf::BarrierFunction, S::StaticData) 
    initialguess = bf.initialguess(S)

    if !bf.isadmissible(initialguess, S)
        throw(ErrorException("Initial guess is not admissible. Please check explicitly provided bounds."))
    end    
    
    return IterationData(initialguess, bf)
end

"""
    set!(I::IterationData, x::AbstractVector{Float64}, S::StaticData)

Sets iteration vector and reassembles dependencies, i.e. gradient and hessian.
"""
function set!(I::IterationData, x::AbstractVector{Float64}, S::StaticData)
    I.x = x
    assemble!(I, S)
end

"""
    apply_descent!(
        I::IterationData,
        t::Union{Float64, Missing},
        c::Union{Vector{Float64}, Missing},
        v::Vector{Float64},
        S::StaticData;
        tracker::DescentTracker = DescentTracker(),
        assemblytracker::AssemblyTracker = AssemblyTracker(),
        backtracking::Bool = true,
        force_nobacktracking::Bool = false
    )

Applies descent direction I.dx to current iterate I.x including potential damping.
Also performs backtracking for long and adaptive schemes
including Armijo line search if required.
Finally, assembles I at new iterate and corresponding derivatives.
Further information is provided in the tracker and asssemblytracker.
"""
function apply_descent!(
    I::IterationData,
    t::Union{Float64, Missing},
    c::Union{Vector{Float64}, Missing},
    v::Vector{Float64},
    S::StaticData;
    tracker::DescentTracker = DescentTracker(),
    assemblytracker::AssemblyTracker = AssemblyTracker(),
    backtracking::Bool = true,
    force_nobacktracking::Bool = false
)
    iterations::Int64 = 0
    r::Float64 = 1.0

    fac::Float64 = 1.0
    if S.damping
        starnormtracker = StarnormTracker()

        lambda = starnorm(v, I.dx)

        if starnormtracker.failed
            fail_damping!(tracker)
            return
        end

        xi = lambda^2 / (1 + lambda)
        fac = 1 / (1 + xi)
    end

    x = I.x - fac * I.dx

    override = ismissing(S.backtracking_override) ? backtracking : S.backtracking_override
    use_backtracking = force_nobacktracking ? false : override

    if use_backtracking
        if S.backtracking_armijo
            val = t * dot(c, I.x) + I.barrier.value(I.x, S)
            grad = t * c .+ I.gradientF
            pred_descent = S.backtracking_armijofactor * fac * grad' * I.dx
        end

        for k = 0:S.backtracking_maxiterations
            if I.barrier.isadmissible(x, S)
                if S.backtracking_armijo
                    valnew = t * dot(c, x) + I.barrier.value(x, S)
                
                    if valnew > val + r * pred_descent
                        r *= S.backtracking_decrementfactor
                    else
                        iterations = k
                        break
                    end
                else
                    iterations = k
                    break
                end
            else
                r *= S.backtracking_decrementfactor
            end
            x = I.x - r * fac * I.dx

            if k == S.backtracking_maxiterations
                fail_direction!(tracker, k, r)
                return
            end
        end
    end

    I.x = x
    assemble!(I, S; tracker=assemblytracker)
    set!(tracker, iterations, r)
end

"""
    starnorm(
        v::AbstractVector{Float64},
        I::IterationData,
        solveLS::Function;
        tracker::StarnormTracker = StarnormTracker()
    )

Computes norm induced by the barrier ``\\Vert v \\Vert^*_x = \\sqrt{v' [F''(x)]^{-1} v}``.
"""
function starnorm(
    v::AbstractVector{Float64},
    I::IterationData,
    solveLS::Function;
    tracker::StarnormTracker = StarnormTracker()
)
    return starnorm(v, solveLS(I.hessianF, v, I.P), tracker=tracker)
end

"""
    starnorm(
        v::AbstractVector{Float64},
        w::AbstractVector{Float64};
        tracker::StarnormTracker = StarnormTracker()
    )

Same as other $(FUNCTIONNAME)(...), but with ``w = [F''(x)]^{-1} v`` precomputed.
"""
function starnorm(
    v::AbstractVector{Float64},
    w::AbstractVector{Float64};
    tracker::StarnormTracker = StarnormTracker()
)
    try 
        return sqrt(v' * w)
    catch e
        if isa(e, DomainError)
            fail!(tracker)
            return NaN
        else
            rethrow()
        end
    end
end

"""
    assemble!(
        I::IterationData,
        S::StaticData;
        tracker::AssemblyTracker = AssemblyTracker()
    )

Assembles iteration values dependent on the present x and static data.
"""
function assemble!(
    I::IterationData,
    S::StaticData;
    tracker::AssemblyTracker = AssemblyTracker()
)
    I.gradientF = I.barrier.gradient(I.x, S)
    hessianF = I.barrier.hessian(I.x,S)

    if tracker.trackhessian
        tracker.hessian = hessianF
    end

    if tracker.trackcondition
        try
            tracker.conditionnumber = cond(hessianF, Inf)
        catch e
            if isa(e, SingularException)
                tracker.conditionnumber = Inf
            else
                rethrow()
            end
        end
    end

    try 
        I.hessianF = S.factorize(hessianF)
    catch e
        if isa(e, PosDefException) || 
            (isa(e, ArgumentError) && occursin("not symmetric", e.msg))

            S.solveLS, S.factorize = select_linearsolver(LU)
            changed_factorization!(tracker)

            try
                I.hessianF = S.factorize(hessianF)
            catch e2 
                if isa(e2, SingularException)
                    I.hessianF = missing
                    hessian_singular!(tracker)
                else
                    rethrow()
                end
            end
        elseif isa(e, SingularException)
            I.hessianF = missing
            hessian_singular!(tracker)
        else
            rethrow()
        end
    end

    try
        I.P = S.computePreconditioner(I.hessianF)
    catch e
        if isa(e, PosDefException)
            S.computePreconditioner = select_preconditioner(LU)
            I.P = S.computePreconditioner(I.hessianF)
            
            changed_factorization!(tracker)
        else
            rethrow()
        end
    end
end
