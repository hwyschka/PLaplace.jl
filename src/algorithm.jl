"""
    solve_plaplace(p::Float64, mesh::Mesh, g::Function, dirichletBoundary::Set{Boundary}; <keyword arguments>) -> PLaplaceData
    
Returns solution to a given problem for the p-Laplace operator with source term and mixed Dirichlet and Neumann boundary conditions. 
Supports either all functions given as continous pointer or discrete evaluations.

# Algorithm

Coming soon. See Sébastien Loisel: Effient algorithms for solving the p-Laplacian in polynomial time

# Arguments
- **Mandatory**
    - `p::Float64`: problem parameter.
    - `mesh::Mesh`: FEM mesh of the domain.
    - `g::Function`: Dirichlet boundary condition.
    - `dirichletBoundary::Set{Boundary}`: part of the boundary to hold Dirichlet condition.
- **Optional**
    - `f::Function`: source term.
    - `neumannBoundary::Set{Boundary}`: part of the boundary to hold Neumann condition.
    - `h::Function`: Neumann boundary condition.
    - `qdim::Int64=1`: image dimension of f, g, h and respectively u, v. (Mandatory for vector-valued solutions.)
    - `eps::Float64=1e-6`: accuracy of main path-following.
    - `stepsize::Stepsize=ADAPTIVE`: stepsize of internal path-following schemes.
    - `maxIter::Int64=1000`: maximal number of iteriations of internal path-following schemes.
    - `kappa::Float64=10.0`: stepsize parameter of long-step path-following variants.
    - `maxIterationsBacktracking::Int64=25`: maximal number of iterations in backtracking of long-step path-following variants.
    - `decrementFactorBacktracking::Float64=0.25`: decrement factor in backtracking of long-step path-following variants.
    - `solver::LinearSolver=CHOLESKY`: solver for linear systems.
    - `preconditioner::Preconditioner=NONE`: preconditioner for linear systems.
    - `useHarmonicProlongation::Bool=true`: prolongation of boundary condition as linear solution or by 0.
    - `objectiveFileName::String=""`: file name for objective log if set.
    - `conditionFileName::String=""`: file name for condition log if set.
    - `consoleOutput::Bool=true`: select if current values are printed to console during the iteration.
"""
function solve_plaplace(p::Float64, mesh::Mesh, g::Function, dirichletBoundary::Set{Boundary}; 
    f::Function=emptyfunction, neumannBoundary::Set{Boundary}=Set{Boundary}(), h::Function=emptyfunction, qdim::Int64=1,
    eps::Float64=1e-6, stepsize::Stepsize=ADAPTIVE, maxIter::Int64=1000, kappa::Float64=10.0, maxIterationsBacktracking::Int64=25, 
    decrementFactorBacktracking::Float64=0.25, solver::LinearSolver=CHOLESKY, preconditioner::Preconditioner=NONE, 
    useHarmonicProlongation::Bool=true, objectiveFileName::String="", conditionFileName::String="", consoleOutput::Bool=true)

    _f = f == emptyfunction ? zeros(Float64, mesh.nnodes * qdim) : evaluate_mesh_function(mesh, f, qdim=qdim) 
    _g = (isempty(dirichletBoundary) || g == emptyfunction) ? zeros(Float64, mesh.nnodes * qdim) : evaluate_mesh_function(mesh, g, dirichletBoundary, qdim=qdim)
    _h = (isempty(neumannBoundary) || h == emptyfunction) ? zeros(Float64, mesh.nnodes * qdim) : evaluate_mesh_function(mesh, h, neumannBoundary, qdim=qdim)
    
    return solve_plaplace(p, mesh, _g, dirichletBoundary, f=_f, neumannBoundary=neumannBoundary, h=_h, 
                    qdim=qdim, eps=eps, stepsize=stepsize,maxIter=maxIter,kappa=kappa, maxIterationsBacktracking=maxIterationsBacktracking, 
                    decrementFactorBacktracking=decrementFactorBacktracking, solver=solver, preconditioner=preconditioner,
                    useHarmonicProlongation=useHarmonicProlongation, objectiveFileName=objectiveFileName, conditionFileName=conditionFileName,
                    consoleOutput=consoleOutput)
end

function solve_plaplace(p::Float64, mesh::Mesh, g::AbstractVector{Float64}, dirichletBoundary::Set{Boundary}; 
    f::AbstractVector{Float64}=Vector{Float64}(), neumannBoundary::Set{Boundary}=Set{Boundary}(), h::AbstractVector{Float64}=Vector{Float64}(), qdim::Int64=1,
    eps::Float64=1e-6, stepsize::Stepsize=ADAPTIVE, maxIter::Int64=1000, kappa::Float64=10.0, maxIterationsBacktracking::Int64=25, 
    decrementFactorBacktracking::Float64=0.25, solver::LinearSolver=CHOLESKY, preconditioner::Preconditioner=NONE, 
    useHarmonicProlongation::Bool=true, objectiveFileName::String="", conditionFileName::String="", consoleOutput::Bool=true)    

    if p < 1 || isinf(p)
        throw(DomainError(p, "This package only supports 1 ≤ p ≤ ∞."))
    end

    times = zeros(Float64, 3)
    iterations = zeros(Int64, 2)

    times[1] += @elapsed outputData = PLaplaceData(mesh, qdim, p, eps, stepsize)
    consoleOutput && print_defaultdata(outputData)

    consoleOutput && println("Initializing static data")
    times[1] += @elapsed staticData = StaticData(p, mesh, dirichletBoundary, neumannBoundary, 
                                                    f, g, h, qdim,
                                                    eps, maxIter, kappa, 
                                                    maxIterationsBacktracking, decrementFactorBacktracking,
                                                    solver, preconditioner, useHarmonicProlongation, 
                                                    consoleOutput, objectiveFileName, conditionFileName)
    add_staticdata!(outputData, staticData)
    consoleOutput && println("Static data assembly finished")

    consoleOutput && println("Initializing iteration data")
    times[1] += @elapsed iterationData = IterationData(compute_initialguess(staticData))
    times[1] += @elapsed assemble!(iterationData, staticData)
    consoleOutput && println("Iteration data assembly finished")

    AuxilliaryPathFollowing, MainPathFollowing = select_pathfollowing(stepsize)
    times[2] += @elapsed iterations[1] += AuxilliaryPathFollowing(staticData, iterationData)
    if iterations[1] > 0
        times[3] += @elapsed iterations[2] += MainPathFollowing(staticData, iterationData)
    end

    add_performance!(outputData, times, iterations)
    add_iterationdata!(outputData, iterationData, staticData)
    consoleOutput && print_statistics(outputData)

    return outputData
end
