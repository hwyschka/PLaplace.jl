"""
$(TYPEDSIGNATURES)

Throws error if given parameters are not in the designated range.
"""
function check_parameters(
    p::Float64,
    qdim::Int64,
    eps::Float64,
    kappa::Float64,
    kappa_updates::Array{Int64},
    kappa_powers::Array{Float64}
)
    msg = ""

    if p < 1 || isinf(p)
        msg *= "\n - This package only supports 1 ≤ p < ∞."
    end

    if qdim < 1
        msg *= "\n - The functions have to have a positive number of components qdim."
    end

    if eps < 0
        msg *= "\n - The algorithm needs positive termination accuracy eps."
    end

    if kappa < eps
        msg *= "\n - The stepsize factor kappa has to be larger than the accuracy eps."
    end

    if length(kappa_updates) != 3
        msg *= "\n - The regions for updates of kappa have to be given by 3 numbers."
    end

    if length(kappa_powers) != 3
        msg *= "\n - 3 powers for the 3 different update regions have to be provided."
    end

    if !isempty(msg)
        errormsg = 
            "\n The submitted problem violates the following parameter restrictions:" * msg
        throw(ErrorException(errormsg))
    end
end

"""
    solve_plaplace(
        p::Float64,
        mesh::Mesh,
        g::Union{AbstractVector{Float64}, Function},
        dirichlet_boundary::Set{Boundary};
        <keyword arguments>
    ) -> PLaplaceData
    
Returns solution to a given problem for the p-Laplace operator with a Dirichlet condition on
(a part of) the boundary. Additionally, a volume source term and Neumann conditions on a
non-intersecting part of the boudary can be specified.
The solution is computed by an interior-point method applied to a reformulated problem.
For more information on the algorithm, please look at the [algorithm](@ref algorithm-theory)
page of the documentation.

Below you find a list of mandatory and keyword arguments indicating which functions can be
specified continously or discrete and which parameters can be adjusted.
Further, the outputs can be controlled in particular if log and vtk files are generated.
It is also possible to track the condition number of the system matrix during the iteration,
but please not that this significantly impacts the performance and requires a lot of memory. 

# Mandatory Arguments
- `p::Float64`: problem parameter.
- `mesh::Mesh`: FEM mesh of the domain.
- `g::Union{AbstractVector{Float64}, Function}`: Dirichlet boundary condition.
- `dirichlet_boundary::Union{Set{Boundary}, Set{Int64}}`: Parts of the physical boundary to
    hold the Dirichlet condition. Either given by set of named boundaries or node ids.

# Keyword Arguments
- `f::Union{AbstractVector{Float64}, Function, Missing} = missing`:
    Discrete or continous source term.
- `neumann_boundary::Union{Set{Boundary}, Set{Int64}, Missing} = missing`:
    Parts of the physical boundary to hold Neumann condition.
    Either given by set of named boundaries or boundary element ids.
- `h::Union{AbstractVector{Float64}, Function, Missing} = missing`:
    Discrete or continous Neumann boundary condition. Requires `neumann_boundary` to be set.
- `boundary_prolongation::Union{AbstractVector{Float64}, Missing} = missing`:
    Prescribed discrete prolongation of the Dirichlet condition to the full domain.
    Will be computed during the algorithm if not set.
- `fixed_nodes::Union{Set{Boundary}, Set{Int64}, Missing} = missing,`:
    Set of fixed node indices, i.e. additional homogeneous Dirichlet nodes.
    Either given by set of named boundaries or node ids.
    Is ignored if `boundary_prolongation` is set.
- `qdim::Int64 = 1`:
    Number of components of the problem, respectively f, g and h as well as u and v.
    Is mandatory to set for vector-valued setting.
- `eps::Float64 = 1e-6`:
    Accuracy of the solution.
    In particular termination criterion for the main path-following.
- `stepsize::Stepsize = ADAPTIVE`:
    Specifier for using an alternative stepping scheme in the internal path-following.
- `maxiterations::Int64 = 1000`:
    Maximal number of iteriations per stage of the internal path-following schemes.
- `kappa::Float64 = 10.0`:
    Parameter used for updates of the step length in long-step path-following variants.
- `backtracking_maxiterations::Int64 = 25`: 
    Maximal number of iterations in backtracking search for the update of the iterate
    in long-step path-following variants.
- `backtracking_decrementfactor::Float64 = 0.25`:
    Decrement factor in backtracking search of long-step path-following variants.
- `solver::LinearSolver = CHOLESKY`:
    Specifier for the solver internally used for linear systems.
- `preconditioner::Preconditioner = NONE`:
    Specifier for the preconditioner internally used for linear systems.
- `useharmonicprolongation::Bool = true`:
    Flag for prolongation of boundary condition as solution to the harmonic problem or by 0.
    Is ignored if `boundary_prolongation` is set.
- `vtkfile::Union{String,Missing} = missing`:
    File name for writing solution to vtk. If not set, no file will be written.
- `verbose::Bool = true`:
    Flag if log is written to console during the iteration.
- `logfile::Union{String,Missing} = missing`:
    File name for writing log to a file. If not set, log will not be exported.
- `logcondition::Bool = false`:
    Flag for tracking condition of the system matrix is tracked and included in log file.
    Will only be exported if `logfile` is set.
"""
function solve_plaplace(
    p::Float64,
    mesh::Mesh,
    g::Union{AbstractVector{Float64}, Function},
    dirichlet_boundary::Union{Set{Boundary}, Set{Int64}}; 
    f::Union{AbstractVector{Float64}, Function, Missing} = missing,
    neumann_boundary::Union{Set{Boundary}, Set{Int64}, Missing} = missing,
    h::Union{AbstractVector{Float64}, Function, Missing} = missing,
    boundary_prolongation::Union{AbstractVector{Float64}, Missing} = missing,
    fixed_nodes::Union{Set{Boundary}, Set{Int64}, Missing} = missing,
    qdim::Int64 = 1,
    eps::Float64 = 1e-6,
    stepsize::Stepsize = ADAPTIVE,
    maxiterations::Int64 = 1000,
    kappa::Float64 = 10.0,
    kappa_updates::Array{Int64} = [2, 8, 15],
    kappa_powers::Array{Float64} = [2.0, 0.5, 0.25],
    backtracking_maxiterations::Int64 = 25, 
    backtracking_decrementfactor::Float64 = 0.25,
    solver::LinearSolver = CHOLESKY,
    preconditioner::Preconditioner = NONE, 
    useharmonicprolongation::Bool = true,
    vtkfile::Union{String,Missing} = missing,
    verbose::Bool = true,
    logfile::Union{String,Missing} = missing,
    logcondition::Bool = false
)
    check_parameters(p, qdim, eps, kappa, kappa_updates, kappa_powers)

    logdata = LogData(verbose, logfile, logcondition)
    log_start(logdata)
    log_defaultdata(logdata, mesh, eps, p, string(stepsize))

    log_setup(logdata, 0)
    _g = g isa Function ? evaluate_mesh_function(mesh, g, dirichlet_boundary, qdim=qdim) : g
    _f = f isa Function ? evaluate_mesh_function(mesh, f, qdim=qdim) : f
    _h = h isa Function ? evaluate_mesh_function(mesh, h, neumann_boundary, qdim=qdim) : h

    _dirichlet_boundary = dirichlet_boundary isa Set{Boundary} ? 
        extract_nodes(dirichlet_boundary) : dirichlet_boundary

    _fixed_nodes = fixed_nodes isa Set{Boundary} ? 
        extract_nodes(fixed_nodes) : fixed_nodes

    _neumann_boundary = neumann_boundary isa Set{Boundary} ? 
        extract_elements(neumann_boundary) : neumann_boundary
    log_setup(logdata, 1)

    algorithmdata = AlgorithmData()

    tsetup::Float64 = 0.0

    tsetup += @elapsed outputdata = PLaplaceData(
        mesh,
        _dirichlet_boundary,
        _fixed_nodes,
        _neumann_boundary,
        _f,
        _h,
        qdim,
        p,
        eps,
        stepsize
    )

    log_setup(logdata, 2)
    tsetup += @elapsed barrierFunction = BarrierFunction()
    log_setup(logdata, 3)

    log_setup(logdata, 4)
    tsetup += @elapsed staticdata = StaticData(
        p,
        mesh,
        _dirichlet_boundary,
        _fixed_nodes,
        _neumann_boundary,
        _f,
        _g,
        _h,
        boundary_prolongation,
        qdim,
        eps,
        maxiterations,
        kappa,
        kappa_updates,
        kappa_powers,
        backtracking_maxiterations,
        backtracking_decrementfactor,
        solver,
        preconditioner,
        useharmonicprolongation
    )
    add_staticdata!(outputdata, staticdata)
    log_setup(logdata, 5)

    log_setup(logdata, 6)
    tsetup += @elapsed iterationdata = IterationData(barrierFunction, staticdata)
    log_setup(logdata, 7)

    algorithmdata.tsetup = tsetup

    AuxilliaryPathFollowing, MainPathFollowing = select_pathfollowing(stepsize)
    algorithmdata.taux = @elapsed AuxilliaryPathFollowing(
        iterationdata,
        algorithmdata,
        staticdata,
        logdata
    )
    if algorithmdata.Naux > 0
        algorithmdata.tmain = @elapsed MainPathFollowing(
            iterationdata, 
            algorithmdata, 
            staticdata,
            logdata
        )
    end

    add_algorithmdata!(outputdata, algorithmdata)
    log_statistics(algorithmdata, logdata)

    if !ismissing(vtkfile)
        write_result_to_vtk(vtkfile, outputdata)
    end

    return outputdata
end
