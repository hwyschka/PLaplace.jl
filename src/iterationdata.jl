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

    "gradient of the barrier function at x"
    gradientF::AbstractVector{Float64}
    "factorized hessian of the barrier function at x"
    hessianF::Any
    "unfactorized hessian (only set if condition needs to be computed)"
    unfactorizedHessianF::SparseMatrixCSC{Float64, Int64}
    "preconditioner for the hessian"
    P::Any

    "message if method fails or changed"
    msg::String
end

function IterationData(x::AbstractVector{Float64})
    le = length(x)
    return IterationData(x, zeros(Float64, le), zeros(Float64, le), spzeros(Float64, le, le), 
                            spzeros(Float64, le, le), LinearAlgebra.I, "-")
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
    apply_descent!(I::IterationData, S::StaticData; useBacktracking=true)

Performs backtracking with descent direction I.dx required for long and adaptive stepsizes
"""
function apply_descent!(I::IterationData, S::StaticData; useBacktracking=true)
    r::Float64 = 1
    x = I.x - r * I.dx

    # Do stepsize backtracking (only for large-stepping schemes)
    maxIter = useBacktracking ? S.maxIterationsBacktracking : 0
    for k = 1:maxIter
        if isadmissible(x,S)
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
    isadmissible(x::AbstractVector{Float64}, S::StaticData) -> Bool
    
Return true if the value x is in the admissible set and false otherwise.
"""
function isadmissible(x::AbstractVector{Float64}, S::StaticData)
    z = zeros(Float64, S.m)
    for (key, val) in S.D
        z -= (val * x[1:S.lengthu] + S.b[key]).^2
    end

    s = x[(S.lengthu + 1):(S.lengthu + S.m)]

    for si in s
        si <= 0 && return false
    end
        
    tau = S.R .- (S.omega .* s)
    for taui in tau
        taui <= 0 && return false
    end

    z += s.^(2 / S.p)
    for zi in z
        zi <= 0 && return false
    end
    
    return true
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
    I.gradientF = compute_gradient_finite(I.x, S)
    hessianF = compute_hessian_finite(I.x,S)

    if S.logCondition
        I.unfactorizedHessianF = hessianF
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
    compute_hessian_finite(x::AbstractVector{Float64}, S::StaticData) -> SparseMatrixCSC{Float64, Int64}

Returns hessian of the barrier function at x for p finite.
"""
function compute_hessian_finite(x::AbstractVector{Float64}, S::StaticData)
    s = x[(S.lengthu + 1):(S.lengthu + S.m)]
    tau = S.R .- (S.omega .* s)

    y = Dict{Tuple{Int64,Int64},AbstractVector{Float64}}()
    for (key, val) in S.D
        y[key] = val * x[1:S.lengthu] + S.b[key]
    end    

    z = zeros(Float64, S.m)
    z += s.^(2 / S.p)
    for (key, val) in y
        z -= val.^2
    end

    F_uu1 = spzeros(Float64,S.lengthu,S.lengthu)
    F_uu2 = spzeros(Float64,S.lengthu,S.lengthu)
    for (key1, val1) in S.D
        F_uu1 += val1' * Diagonal(z.^(-1)) * val1
        for (key2, val2) in S.D
            F_uu2 += (Diagonal(y[key1]) * val1)' * Diagonal(z.^(-2)) * (Diagonal(y[key2]) * val2)
        end
    end
    F_uu = 2 * F_uu1 + 4 * F_uu2
    
    Fus = spzeros(S.lengthu, S.m)
    for (key, val) in S.D
        Fus -= (Diagonal(y[key]) * val)' * Diagonal(z.^(-2)) * Diagonal(s.^(2.0/S.p - 1))
    end
    Fus *= 4 / S.p

    F_ss = zeros(S.m)
    F_ss -= 2 / S.p * (2.0 / S.p - 1.0) .* z.^(-1) .* s.^(2.0 / S.p - 2)
    F_ss += 4 / S.p^2 .* z.^(-2) .* s.^(4.0 / S.p - 2)
    F_ss += S.sigma * s.^(-2)
    F_ss += S.omega.^2 .* tau.^(-2)
    
    hessF = [Symmetric(F_uu) Fus; Fus' Diagonal(F_ss)]

    return hessF
end

function compute_gradient_finite(x::AbstractVector{Float64}, S::StaticData)
    s = x[(S.lengthu + 1):(S.lengthu + S.m)]
    tau = S.R .- (S.omega .* s)

    y = Dict{Tuple{Int64,Int64},AbstractVector{Float64}}()
    for (key, val) in S.D
        y[key] = val * x[1:S.lengthu] + S.b[key]
    end

    z = zeros(Float64, S.m)
    z += s.^(2 / S.p)
    for (key, val) in y
        z -= val.^2
    end

    Fu = zeros(Float64, S.lengthu)
    for (key, val) in S.D
        Fu += val' * (y[key] ./ z)
    end
    Fu *= 2

    Fs = zeros(Float64, S.m)
    Fs -= 2 / S.p * z.^(-1) .* s.^(2.0 / S.p - 1) 
    Fs += S.omega ./ tau 
    Fs -= S.sigma ./ s

    return [Fu; Fs]
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
