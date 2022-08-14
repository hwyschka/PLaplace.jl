"""
    mutable struct StaticData

Stores all data that is static throughout iteration of the interior-point method.
"""
mutable struct StaticData  
    "PDE parameter"
    p::Float64
    "domain dimension"
    d::Int8
    "image dimension"
    qdim::Int64
    "number of nodes"
    n::Int64
    "number of elements"
    m::Int64

    "set of actually calculated nodes"
    calculationNodes::Set{Int64}
    "number of calculated nodes"
    nCalc::Int64
    "number of calculated coefficients"
    lengthu::Int64

    "system vector for convex optimization problem"
    c::AbstractVector{Float64}
    "upper bound for auxilliary variable"
    R::Float64

    "vector of element weights"
    omega::AbstractVector{Float64}
    "discrete derivative matrix"
    D::AbstractDict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}
    "prolonged boundary"
    g::AbstractVector{Float64}
    "discrete derivative of the prolonged boundary condition"
    b::AbstractDict{Tuple{Int64,Int64},AbstractVector{Float64}}

    "constant for barrier parameter"
    sigma::Int8
    "algorithm constant given by Nesterov as 1/9"
    beta::Float64
    "algorithm constant given by Nesterov as 5/36"
    gamma::Float64
    "accuracy"
    eps::Float64
    "inverse tolerance"
    tolInv::Float64
    "accuracy to tolerance factor"
    tolFactor::Float64
    "maximal number of iterations in path-following schemes"
    maxIter::Int64
    "initial stepsize for long and adaptive path-following"
    kappa::Float64
    "maximal number of iterations for backtracking in long and adaptive path-following"
    maxIterationsBacktracking::Int64
    "decrement factor for backtracking in long and adaptive path-following"
    decrementFactorBacktracking::Float64
    
    "function pointer on solving routine"
    solveLS::Function
    "funciton pointer on factorization routine"
    factorize::Function
    "function pointer on preconditioner"
    computePreconditioner::Function
    
    "function pointer on console output routine"
    log::Function
    "switch for logging the objective to file during the iteration"
    logObjective::Bool
    "switch for logging the condition of the linear system to file during the iteration"
    logCondition::Bool
    "file name for objective logging during the auxilliary path-following"
    auxObjectiveFileName::String
    "file name for objective logging during the main path-following"
    mainObjectiveFileName::String
    "file name for condition logging during the auxilliary path-following"
    auxConditionFileName::String
    "file name for condition logging during the main path-following"
    mainConditionFileName::String
end

function StaticData(p::Float64, mesh::Mesh, dirichletBoundary::Set{Boundary}, neumannBoundary::Set{Boundary}, 
                    f::AbstractVector{Float64}, g::AbstractVector{Float64}, h::AbstractVector{Float64}, qdim::Int64,
                    eps::Float64, maxIter::Int64, kappa::Float64, maxIterationsBacktracking::Int64, decrementFactorBacktracking::Float64,
                    solver::LinearSolver, preconditioner::Preconditioner, useHarmonicProlongation::Bool, consoleOutput::Bool,
                    objectiveFileName::String, conditionFileName::String)

    log_console = consoleOutput ? println : emptyfunction
    logObjective = !isempty(objectiveFileName)
    logCondition = !isempty(conditionFileName)
    auxObjectiveFileName = objectiveFileName * "_aux"
    mainObjectiveFileName = objectiveFileName * "_main"
    auxConditionFileName = conditionFileName * "_aux"
    mainConditionFileName = conditionFileName * "_main"
    log_console("   Log files set")

    solveLS, factorize = select_linearsolver(solver)
    computePreconditioner = select_preconditioner(preconditioner)
    log_console("   Linear solver set")

    sigma = compute_sigma(p)
    nu = mesh.nelems * (sigma + 2)
    beta = 1/9
    gamma = 5/36
    tolFactor = (nu + (((beta + sqrt(nu)) * beta) / (1 - beta))) 
    tolInv =  tolFactor / eps
    log_console("   Algorithm parameters computed")

    dirichletNodes = extract_nodes(dirichletBoundary)
    calculationNodes = setdiff(Set(1:mesh.nnodes), dirichletNodes)

    neumannElements = Set{Int64}()
    neumannNodes = Set{Int64}()
    hasNeumannBoundary::Bool = false
    if(!isempty(neumannBoundary))
        hasNeumannBoundary = true
        neumannElements = extract_elements(neumannBoundary)
        neumannNodes = extract_nodes(neumannBoundary)
    end
    hasSourceTerm = any(o -> o != 0, f)
    log_console("   Boundary sets computed")

    if useHarmonicProlongation
        gProl = compute_prolongation_harmonic(g, mesh, hasSourceTerm, f, dirichletNodes, hasNeumannBoundary, h, neumannElements, solveLS, computePreconditioner, qdim=qdim)
    else
        gProl = compute_prolongation_zero(g, mesh)
    end
    log_console("   Boundary prolongation finished")

    omega = assemble_weightmultivector(mesh, qdim=1, order=1)
    D = assemble_derivativetensor(mesh, qdim=qdim)
    Dmod = assemble_derviativetensor_modified(D, dirichletNodes, qdim=qdim)
    b = compute_derivative(D, gProl)
    log_console("   Discrete Derivatives assembled")

    c = compute_systemvector(mesh, p, omega, dirichletNodes, hasNeumannBoundary, neumannElements, h, hasSourceTerm, f, qdim)
    R = compute_upperbound(mesh, p, b, omega, hasNeumannBoundary, neumannElements, h, hasSourceTerm, f, qdim)
    log_console("   System values set")

    return StaticData(p, mesh.d, qdim, mesh.nnodes, mesh.nelems, 
                        calculationNodes, length(calculationNodes), length(calculationNodes)*qdim,
                        c, R,
                        omega, Dmod, gProl, b,
                        sigma, beta, gamma, eps, tolInv, tolFactor, maxIter, 
                        kappa, maxIterationsBacktracking, decrementFactorBacktracking,
                        solveLS, factorize, computePreconditioner,
                        log_console, logObjective, logCondition,
                        auxObjectiveFileName, mainObjectiveFileName, auxConditionFileName, mainConditionFileName)
end

"""
    compute_sigma(p::Float64) -> Int64
    
Returns factor σ for the barrier function parameter depending on p.
"""
function compute_sigma(p::Float64)
    return p >= 2 ? 1 : 2
end

"""
    compute_systemvector(mesh::Mesh, p::Float64, omega::AbstractVector{Float64}, dirichletNodes::Set{Int64}, 
                            hasNeumannBoundary::Bool, neumannElements::Set{Int64}, h::AbstractVector{Float64},
                            hasSourceTerm::Bool, f::AbstractVector{Float64}, qdim::Int64) -> Vector{Float64}
    
Computes system vector c for a p-Poisson problem to solve with a barrier method.
"""
function compute_systemvector(mesh::Mesh, p::Float64, omega::AbstractVector{Float64}, dirichletNodes::Set{Int64}, 
                                hasNeumannBoundary::Bool, neumannElements::Set{Int64}, h::AbstractVector{Float64},
                                hasSourceTerm::Bool, f::AbstractVector{Float64}, qdim::Int64)
    
    elementsToDrop = Set{Int64}()
    if qdim == 1
        elementsToDrop = dirichletNodes
    else
        for node in dirichletNodes
            union!(elementsToDrop, qdim*(node-1)+1:qdim*node)
        end
    end

    rhs = zeros(Float64, qdim*(mesh.nnodes-length(dirichletNodes)))
    
    if hasSourceTerm
        if length(f) == mesh.nnodes*qdim
            M = assemble_massmatrix(mesh, qdim=qdim, order=3)
            rhs -= (M * f)[1:end .∉ [elementsToDrop]]
        elseif mod(length(f), mesh.nelems*qdim) == 0
            nPoints = length(f) / (mesh.nelems * qdim)
            quadOrder = quadrature_order(mesh.d, nPoints)
            
            E = assemble_basismatrix(mesh, qdim=qdim, order=quadOrder)
            W = Diagonal(assemble_weightmultivector(mesh, qdim=qdim, order=quadOrder))
            rhs -= (E' * W * f)[1:end .∉ [elementsToDrop]]
        else
            throw(DomainError(f,"Dimension Missmatch"))
        end
    end
    
    if hasNeumannBoundary
        if length(h) == mesh.nnodes*qdim
            N = assemble_massmatrix_boundary(mesh, boundaryElements=neumannElements, qdim=qdim, order=3)
            rhs -= (N * h)[1:end .∉ [elementsToDrop]]
        elseif mod(length(h),mesh.nboundelems*qdim) == 0
            nPoints = div(length(h), mesh.nboundelems * qdim)
            quadOrder = quadrature_order(mesh.d-1, nPoints)

            E = assemble_basismatrix_boundary(mesh, boundaryElements=neumannElements, qdim=qdim, order=quadOrder)
            W = Diagonal(assemble_weightmultivector_boundary(mesh, qdim=qdim, order=quadOrder))
            rhs -= (E' * W * h)[1:end .∉ [elementsToDrop]]
        else
            throw(DomainError(h,"Dimension Missmatch"))
        end
    end

    return [rhs; omega ./ p]
end

"""
    compute_upperbound(mesh::Mesh, p::Float64, b::AbstractDict{Tuple{Int64,Int64},AbstractVector{Float64}}, omega::AbstractVector{Float64}, 
                        hasNeumannBoundary::Bool, neumannElements::Set{Int64}, h::AbstractVector{Float64}, 
                        hasSourceTerm::Bool, f::AbstractVector{Float64}, qdim::Int64) -> Float64

Computes constant upper bound R on the element derivative.
"""
function compute_upperbound(mesh::Mesh, p::Float64, b::AbstractDict{Tuple{Int64,Int64},AbstractVector{Float64}}, omega::AbstractVector{Float64}, 
                            hasNeumannBoundary::Bool, neumannElements::Set{Int64}, h::AbstractVector{Float64}, 
                            hasSourceTerm::Bool, f::AbstractVector{Float64}, qdim::Int64)
    L::Float64 = stripwidth(mesh)
    
    normg::Float64 = xpnorm(p, b, omega)
    normf::Float64 = 0
    normh::Float64 = 0
    
    if hasSourceTerm
        normf = qnorm(p, f, mesh, qdim=qdim, order=5)
    end

    if hasNeumannBoundary
        normh = qnorm_boundary(p, h, mesh, boundaryElements=neumannElements, qdim=qdim, order=5)
    end

    if p == 1
        R = 2 + 2 * normg / (1 - L * normf)
    else
        q = conjugated_exponent(p)
        R = 2 + 8 * normg^(p) + 4 * L^q * (p / 2)^(1 / (1 - p)) * (p - 1) * normf
    end
    
    return R
end

"""
    compute_prolongation_harmonic(g::Function, mesh::Mesh, hasSourceTerm::Bool, f::Function, 
    dirichletNodes::Set{Int64}, hasNeumannBoundary::Bool, h::Function, neumannElements::Set{Int64}, neumannNodes::Set{Int64}, 
    solveLS::Function, preconditioner::Function; qdim=1) -> Vector{Float64}
    
Returns solution of linear Laplace problem to prolongate boundary condition g to the domain.

# Arguments
- `g::AbstractVector{Float64}`: problem parameter.
- `mesh::Mesh`: mesh of the domain.
- `hasSourceTerm::Bool`: switch for evaluating source term.
- `f::AbstractVector{Float64}`: source term.
- `dirichletNodes::Set{Int64}`: nodes of the boundary to hold Dirichlet condition.
- `hasNeumannBoundary::Bool`: switch for evaluating Neumann boundary.
- `h::AbstractVector{Float64}`: Neumann boundary condition.
- `neumannElements::Set{Int64}`: edges of the boundary to hold Neumann condition.
- `solveLS::Function`: function pointer to solve linear system.
- `preconditioner::Function`: function pointer to compute preconditioner for linear system.
- `qdim::Int64=1`: image dimension of f, h and g.

"""
function compute_prolongation_harmonic(g::AbstractVector{Float64}, mesh::Mesh, hasSourceTerm::Bool, f::AbstractVector{Float64}, 
                                        dirichletNodes::Set{Int64}, hasNeumannBoundary::Bool, h::AbstractVector{Float64}, neumannElements::Set{Int64}, 
                                        solveLS::Function, preconditioner::Function; qdim::Int64=1)

    rhs = zeros(Float64, qdim*mesh.nnodes)

    if hasSourceTerm
        if length(f) == mesh.nnodes*qdim
            M = assemble_massmatrix(mesh, qdim=qdim, order=3)
            rhs += M * f
        elseif mod(length(f), mesh.nelems*qdim) == 0
            nPoints = length(f) / (mesh.nelems * qdim)
            quadOrder = quadrature_order(mesh.d, nPoints)

            E = assemble_basismatrix(mesh, qdim=qdim, order=quadOrder)
            W = Diagonal(assemble_weightmultivector(mesh, qdim=qdim, order=quadOrder))
            rhs += E' * W * f
        else
            throw(DomainError(f,"Dimension Missmatch"))
        end
    end
    
    if hasNeumannBoundary
        if length(h) == mesh.nnodes*qdim
            N = assemble_massmatrix_boundary(mesh, boundaryElements=neumannElements, qdim=qdim, order=3)
            rhs += N * h
        elseif mod(length(h), mesh.nboundelems*qdim) == 0
            nPoints = div(length(h), mesh.nboundelems * qdim)
            quadOrder = quadrature_order(mesh.d-1, nPoints)

            E = assemble_basismatrix_boundary(mesh, boundaryElements=neumannElements, qdim=qdim, order=quadOrder)
            W = Diagonal(assemble_weightmultivector_boundary(mesh, qdim=qdim, order=quadOrder))
            rhs += E' * W * h
        else
            throw(DomainError(h,"Dimension Missmatch"))
        end
    end
    
    L = assemble_laplacian(mesh, qdim=qdim)
    assemble_dirichletcondition!(L, dirichletNodes, rhs=rhs, bc=g, qdim=qdim)

    return solveLS(L, rhs, preconditioner(L))
end

"""
    compute_prolongation_zero(g::AbstractVector{Float64}, mesh::Mesh; qdim::Int64=1) -> Vector{Float64}
    
Returns discrete prolongation of g to the whole domain by zero on every node.
"""
function compute_prolongation_zero(g::AbstractVector{Float64}, mesh::Mesh; qdim::Int64=1)
    return zeros(Float64, qdim*mesh.nnodes) .+ g
end

"""
    compute_initialguess(S::StaticData) -> Vector{Float64}
    
Returns start vector for an auxilliary path-following scheme.
"""
function compute_initialguess(S::StaticData)
    t = zeros(Float64, S.m)
    for (key, val) in S.b
        t += val.^2
    end

    s = 1 .+ t.^(S.p/2)

    return [zeros(Float64, S.lengthu); s]
end

"""
    add_staticdata!(data::PLaplaceData, S::StaticData)

Adds required values from the static data to the output data, especially the result.
"""
function add_staticdata!(data::PLaplaceData, S::StaticData)
    data.g = S.g
    data.calculationNodes = S.calculationNodes
end
