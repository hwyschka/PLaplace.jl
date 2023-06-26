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
    "condition of the hessian (only set if condition needs to be computed)"
    conditionHessian::Union{Missing,Float64}
    "preconditioner for the hessian"
    P::Any

    "message if method fails or changed"
    msg::String
end

function IterationData(x::AbstractVector{Float64}, bf::BarrierFunction)
    le = length(x)
    IterationData(
        x,
        zeros(Float64, le),
        bf,
        zeros(Float64, le),
        spzeros(Float64, le, le), 
        missing,
        LinearAlgebra.I,
        "-"
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
function apply_descent!(I::IterationData, S::StaticData; useBacktracking=true)
    r::Float64 = 1
    x = I.x - r * I.dx

    # Do stepsize backtracking (only for large-stepping schemes)
    maxIter = useBacktracking ? S.maxIterationsBacktracking : 0
    for k = 1:maxIter
        if I.barrier.isadmissible(x,S)
            S.log("r = $r in domain ")
            break
        else
            S.log("r = $r not in domain")
            r *= S.decrementFactorBacktracking
        end
        x = I.x - r * I.dx
    end

    I.x = x
    assemble!(I, S)
end

"""
    starnorm(v::AbstractVector{Float64}, I::IterationData, solveLS::Function)
    starnorm(v::AbstractVector{Float64}, w::AbstractVector{Float64})

Computes ||v||^*_x = sqrt(v' F''(x)^(-1) v) potentially with w = F''(x)^(-1) v precomputed.
"""
function starnorm(v::AbstractVector{Float64}, I::IterationData, solveLS::Function)
    return sqrt(v' * solveLS(I.hessianF, v, I.P))
end

function starnorm(v::AbstractVector{Float64}, w::AbstractVector{Float64})
    return sqrt(v' * w)
end

"""
    assemble!(I::IterationData, S::StaticData)

Assembles iteration values dependend on the present x and static data.
"""
function assemble!(I::IterationData, S::StaticData)
    I.gradientF = I.barrier.gradient(I.x, S)
    hessianF = I.barrier.hessian(I.x,S)

    if S.logCondition
        I.conditionHessian = cond(hessianF, Inf)
    end

    try 
        I.hessianF = S.factorize(hessianF)
    catch e
        if isa(e, PosDefException)
            S.solveLS, S.factorize = select_linearsolver(LU)
            I.hessianF = S.factorize(hessianF)

            I.msg *= " Changed to LU fact. - "
            S.log("Hessian not numerical pos def.\n Changed factorization to LU.")
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

            I.msg *= " Changed to LU prec. - "
            S.log("Hessian not numerical pos def.\n Changed precondtioner to LU.")
        else
            rethrow()
        end
    end
end

"""
    add_iterationdata!(data::PLaplaceData, I::IterationData, S::StaticData)

Adds required values from the iteration data to the output data, especially the result.
"""
function add_iterationdata!(data::PLaplaceData, I::IterationData, S::StaticData)
    data.u = I.x[1:S.lengthu]
    assemble_result!(data)
    data.msg = I.msg
    data.eps = S.eps
end
