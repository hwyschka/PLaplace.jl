"""
$(TYPEDEF)

Public type used to specify the linear solver in the solution algorithm.

# Available Options
- `LU`:
    Direct solver with LU decomposition.
- `CHOLESKY`:
    Direct solver with Cholesky decomposition.
    Default options since Hessian is supposed to be symmetric positive definite.
- `CG`:
    Iterative solver using the conjugate gradients methods.
- `BICGSTAB`:
    Iterative solver using the bi-conjugate gradients stabilized method.
- `GMRES`:
    Iterative solver using the generalized minimal residual method.
    Has proven convergence for unsymmetric systems, but requires a lot of memory.
"""
@enum LinearSolver begin
    LU
    CHOLESKY
    CG
    BICGSTAB
    GMRES
end

"""
$(TYPEDEF)

Public type used to specify a preconditioner for the linear systems in the solution process.

# Available Options
- `NONE`:
    Dummy option for no preconditioner. In particular used with direct solvers.
- `JACOBI`:
    Preconditioner with diagonal matrix.
- `ILU`:
    Preconditioner with incomplete LU decomposition.
- `ICHOLESKY`:
    Preconditioner with incomplete Cholesky decomposition.
- `AMG`:
    Preconditioner with algebraic multi-grid.
"""
@enum Preconditioner begin
    NONE
    JACOBI
    ILU
    ICHOLESKY
    AMG
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for solving linear system with julias default \\ operator.
Usually intended to be called not with a matrix A itself
but a factorization, e.g. LU or Cholesky.
"""
function solve_default(
    A::Any,
    b::AbstractVector{Float64},
    P::Any
)
    return A\b
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for solving linear system with a preconditioned conjugated gradients method.
"""
function solve_cg(
    A::SparseMatrixCSC{Float64,Int64},
    b::AbstractVector{Float64},
    P::Any
)
    return cg(A, b, Pl=P)
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for solving linear system with a preconditioned stabilized
bi-conjugated gradients method.
"""
function solve_bicgstab(
    A::SparseMatrixCSC{Float64,Int64},
    b::AbstractVector{Float64},
    P::Any
)
    return bicgstabl(A, b, Pl=P)
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for solving linear system with a preconditioned GMRES.
"""
function solve_gmres(
    A::SparseMatrixCSC{Float64,Int64},
    b::AbstractVector{Float64},
    P::Any
)
    return gmres(A,b,Pl=P)
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for computing a preconditioner for a matrix A or its factorization without effect,
i.e. returns identity.
"""
function preconditioner_none(A::Any)
    return LinearAlgebra.I
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for computing a Jacobi (/Diagonal) preconditioner for a matrix A.
"""
function preconditioner_jacobi(A::SparseMatrixCSC{Float64,Int64})
    return DiagonalPreconditioner(A)
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for computing an incomplete LU preconditioner for a matrix A.
"""
function preconditioner_ilu(A::SparseMatrixCSC{Float64,Int64})
    return ilu(A, τ = 0.1)
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for computing an incomplete Cholesky preconditioner for a matrix A.
"""
function preconditioner_illt(A::AbstractSparseMatrix{Float64})
    return CholeskyPreconditioner(A)
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for computing an algebraic multi-grid preconditioner for a matrix A.
"""
function preconditioner_amg(A::AbstractSparseMatrix{Float64})
    return AMGPreconditioner(A)
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for not computing a factorization of a matrix A, but keeping the matrix itself.
Can be used as do nothing factorization for iterative linear solvers.
"""
function factorization_none(A::SparseMatrixCSC{Float64,Int64})
    return A
end

LinearSolvers = Dict{LinearSolver, Function}(LU => solve_default,
                                                CHOLESKY => solve_default,
                                                CG => solve_cg,
                                                BICGSTAB => solve_bicgstab,
                                                GMRES => solve_gmres)

Factorizations = Dict{LinearSolver, Function}(LU => lu,
                                                CHOLESKY => cholesky,
                                                CG => factorization_none,
                                                BICGSTAB => factorization_none,
                                                GMRES => factorization_none)

Preconditioners = Dict{Preconditioner,Function}(NONE => preconditioner_none,
                                                JACOBI => preconditioner_jacobi,
                                                ILU => preconditioner_ilu,
                                                ICHOLESKY => preconditioner_illt,
                                                AMG => preconditioner_amg)

"""
$(TYPEDSIGNATURES)
    
Returns pair of functions of the actual solving routine
and the factorization routine for the given linear solver.
"""
function select_linearsolver(solver::LinearSolver) :: Tuple{Function, Function}
    return get(LinearSolvers, solver, solve_default),
        get(Factorizations, solver, factorization_none)
end

"""
$(TYPEDSIGNATURES)
    
Returns function for computing the given preconditioner.
"""
function select_preconditioner(preconditioner::Preconditioner) :: Function
    return get(Preconditioners, preconditioner, preconditioner_none)
end
