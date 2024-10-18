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

IterationData(bf::BarrierFunction, S::StaticData) = IterationData(bf.initialguess(S), bf)

"""
    set!(I::IterationData, x::AbstractVector{Float64}, S::StaticData)

Sets iteration vector and reassembles dependencies, i.e. gradient and hessian.
"""
function set!(I::IterationData, x::AbstractVector{Float64}, S::StaticData)
    I.x = x
    assemble!(I, S)
end

"""
    apply_descent!(I::IterationData, S::StaticData; useBacktracking=true)

Performs backtracking with descent direction I.dx required for long and adaptive stepsizes
"""
function apply_descent!(
    I::IterationData,
    S::StaticData;
    tracker::DescentTracker = DescentTracker(0,0.0),
    assemblytracker::AssemblyTracker = AssemblyTracker(),
    useBacktracking::Bool = true
)
    iterations::Int64 = 0
    r::Float64 = 1
    x = I.x - r * I.dx

    maxIter = useBacktracking ? S.maxIterationsBacktracking : -1
    for k = 0:maxIter
        if I.barrier.isadmissible(x,S)
            iterations = k
            break
        else
            r *= S.decrementFactorBacktracking
        end
        x = I.x - r * I.dx
    end

    I.x = x
    assemble!(I, S; tracker=assemblytracker)
    set!(tracker, iterations, r)
end

"""
$(TYPEDSIGNATURES)

Computes norm induced by the barrier ``\\Vert v \\Vert^*_x = \\sqrt{v' [F''(x)]^{-1} v}``.
"""
function starnorm(v::AbstractVector{Float64}, I::IterationData, solveLS::Function)
    return sqrt(v' * solveLS(I.hessianF, v, I.P))
end

"""
$(TYPEDSIGNATURES)

Same as other $(FUNCTIONNAME)(...), but with ``w = [F''(x)]^{-1} v`` precomputed.
"""
function starnorm(v::AbstractVector{Float64}, w::AbstractVector{Float64})
    return sqrt(v' * w)
end

"""
    assemble!(I::IterationData, S::StaticData)

Assembles iteration values dependend on the present x and static data.
"""
function assemble!(
    I::IterationData,
    S::StaticData;
    tracker::AssemblyTracker=AssemblyTracker()
)
    I.gradientF = I.barrier.gradient(I.x, S)
    hessianF = I.barrier.hessian(I.x,S)

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
