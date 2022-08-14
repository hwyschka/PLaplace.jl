"""
    @enum LinearSolver

Variable used to specify the linear solver in the solution algorithm.
"""
@enum LinearSolver begin
    LU
    CHOLESKY
    CG
    GMRES
end

"""
    @enum Preconditioner

Variable used to specify a preconditioner for the linear systems in the solution process.
"""
@enum Preconditioner begin
    NONE
    JACOBI
    ILU
    ICHOLESKY
    AMG
end

"""
    solve_default(A::Any, b::AbstractVector, P::Any) -> Vector{Float64}
    
Wrapper for solving linear system with julias default \\ operator.
Usually intended to be called not with a matrix A itself but a factorization, e.g. LU or Cholesky.
"""
function solve_default(A::Any, b::AbstractVector, P::Any)
    return A\b
end

"""
    solve_cg(A::AbstractSparseMatrix{Float64}, b::AbstractVector, P::Any) -> Vector{Float64}
    
Wrapper for solving linear system with a preconditioned conjugated gradients method.
"""
function solve_cg(A::AbstractSparseMatrix{Float64}, b::AbstractVector, P::Any)
    return cg(A,b,Pl=P)
end

"""
    solve_gmres(A::AbstractSparseMatrix{Float64}, b::AbstractVector, P::Any) -> Vector{Float64}
    
Wrapper for solving linear system with a preconditioned GMRES.
"""
function solve_gmres(A::AbstractSparseMatrix{Float64}, b::AbstractVector, P::Any)
    return gmres(A,b,Pl=P)
end

"""
    preconditioner_none(A::Any) -> LinearAlgebra.UniformScaling
    
Wrapper for computing a preconditioner for a matrix A or its factorization without effect, i.e. returns identity.
"""
function preconditioner_none(A::Any)
    return LinearAlgebra.I
end

"""
    preconditioner_jacobi(A::AbstractSparseMatrix{Float64}) -> AbstractPreconditioner
    
Wrapper for computing a Jacobi (/Diagonal) preconditioner for a matrix A.
"""
function preconditioner_jacobi(A::AbstractSparseMatrix{Float64})
    return DiagonalPreconditioner(A)
end

"""
    preconditioner_ilu(A::AbstractSparseMatrix{Float64})
    
Wrapper for computing an incomplete LU preconditioner for a matrix A.
"""
function preconditioner_ilu(A::AbstractSparseMatrix{Float64})
    return ilu(A, τ = 0.1)
end

"""
    ICHOLPrecond(A::AbstractSparseMatrix{Float64})
    
Wrapper for computing an incomplete Cholesky preconditioner for a matrix A.
"""
function preconditioner_illt(A::AbstractSparseMatrix{Float64})
    return CholeskyPreconditioner(A)
end

"""
    preconditioner_amg(A::AbstractSparseMatrix{Float64})
    
Wrapper for computing an algebraic multi-grid preconditioner for a matrix A.
"""
function preconditioner_amg(A::AbstractSparseMatrix{Float64})
    return AMGPreconditioner(A)
end

"""
    factorization_none(A::AbstractSparseMatrix{Float64}) -> AbstractSparseMatrix{Float64}
    
Wrapper for not computing a factorization of a matrix A, but keeping the matrix itself.
"""
function factorization_none(A::AbstractSparseMatrix{Float64})
    return A
end

LinearSolvers = Dict{LinearSolver, Function}(LU => solve_default,
                                                CHOLESKY => solve_default,
                                                CG => solve_cg,
                                                GMRES => solve_gmres)

Factorizations = Dict{LinearSolver, Function}(LU => lu,
                                                CHOLESKY => cholesky,
                                                CG => factorization_none,
                                                GMRES => factorization_none)

Preconditioners = Dict{Preconditioner,Function}(NONE => preconditioner_none,
                                                JACOBI => preconditioner_jacobi,
                                                ILU => preconditioner_ilu,
                                                ICHOLESKY => preconditioner_illt,
                                                AMG => preconditioner_amg)

"""
    select_linearsolver(solver::LinearSolver) -> Function, Function
    
Returns pair of functions of the actual solving routine and the factorization routine for the given linear solver.
"""
function select_linearsolver(solver::LinearSolver)
    return get(LinearSolvers, solver, solve_default), get(Factorizations, solver, factorization_none)
end

"""
    select_preconditioner(preconditioner::Preconditioner) -> Function
    
Returns function for computing the given preconditioner.
"""
function select_preconditioner(preconditioner::Preconditioner)
    return get(Preconditioners, preconditioner, preconditioner_none)
end
