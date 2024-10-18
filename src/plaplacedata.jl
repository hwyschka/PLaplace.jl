"""
$(TYPEDEF)

Main output type for solution of a p-Laplace problem.
Stores all parameters and results obtained by the solution of a pLaplace problem.
Can be used for [Statistics](@ref), [Error Analysis](@ref)
and writing the solution to a .vtk-file.

# Fields
$(TYPEDFIELDS)
"""
mutable struct PLaplaceData
    "Mesh of computational domain."
    mesh::Mesh

    "Dirichlet boundary of the computational domain."
    dirichlet_nodes::Set{Int64}

    "Fixed nodes on the boundary of the computational domain."
    fixed_nodes::Union{Set{Int64}, Missing}

    "Neumann boundary of the computational domain."
    neumann_elements::Union{Set{Int64}, Missing}

    "Domain source term."
    f::Union{AbstractVector{Float64}, Missing}

    "Neumann boundary condition."
    h::Union{AbstractVector{Float64}, Missing}
    
    "Number of components of the solution."
    qdim::Int64

    "PDE parameter."
    p::Float64

    "Upper bound for the gradient of the solution."
    gradient_bound::Union{Float64, Missing}

    "Accuracy of the solution in a variational sense."
    eps::Float64

    "Stepping scheme of the interior point method."
    stepsize::Stepsize
    
    "Number of iterations in the auxilliary path-following.
        Is missing if algorithm did not reach auxilliary stage."
    Naux::Union{Int64,Missing}

    "Number of iterations for the main path-following
        Is missing if algorithm did not reach main stage."
    Nmain::Union{Int64,Missing}

    "Time required for the setup
        Is missing if algorithm did not reach setup stage."
    tsetup::Union{Float64,Missing}

    "Time required for the auxilliary path-following
        Is missing if algorithm did not reach auxilliary stage."
    taux::Union{Float64,Missing}

    "Time required for the main path-following
        Is missing if algorithm did not reach main stage."
    tmain::Union{Float64,Missing}

    "Nodes on which the solution was computed, i.e. nodes without Dirichlet condtion.
        Is missing if algorithm failed during setup stage."
    calculation_nodes::Union{Set{Int64},Missing}

    "Dirichlet boundary condition prolonged to the whole domain.
        Is missing if algorithm failed during setup stage."
    g::Union{AbstractVector{Float64},Missing}

    "Shifted solution on calculation nodes without prolonged boundary condition.
        Is missing if algorithm did not converge."
    u::Union{AbstractVector{Float64},Missing}

    "Solution to the p-Laplace problem.
            Is missing if algorithm did not converge."
    v::Union{AbstractVector{Float64},Missing}

    "Notifications from the iteration. In particular contains changes of solvers and
        preconditioners and information on early stops."
    msg::String
end

PLaplaceData(
    mesh::Mesh,
    dirichlet_boundary::Set{Int64},
    fixed_nodes::Union{Set{Int64}, Missing},
    neumann_boundary::Union{Set{Int64}, Missing},
    f::Union{AbstractVector{Float64}, Missing},
    h::Union{AbstractVector{Float64}, Missing},
    qdim::Int64,
    p::Float64,
    eps::Float64,
    stepsize::Stepsize
) = PLaplaceData(
    mesh,
    dirichlet_boundary,
    fixed_nodes,
    neumann_boundary,
    f,
    h,
    qdim,
    p,
    missing,
    eps,
    stepsize,
    missing,
    missing,
    missing,
    missing,
    missing,
    missing,
    missing,
    missing,
    missing,
    "-"
)

"""
$(TYPEDSIGNATURES)

Assembles combined result vector v = u + g for truncated u vector u_I.
"""
function assemble_result!(data::PLaplaceData)
    if ismissing(data.u)
        return
    end

    n = data.mesh.nnodes
    v = zeros(data.qdim*n)
    k = 1
    for i = 1:n
        v[data.qdim*(i-1)+1:data.qdim*i] = data.g[data.qdim*(i-1)+1:data.qdim*i]
        
        if i in data.calculation_nodes
            v[data.qdim*(i-1)+1:data.qdim*i] += data.u[data.qdim*(k-1)+1:data.qdim*k]
            k += 1
        end
    end
    
    data.v = v
end

"""
$(TYPEDSIGNATURES)

Returns if the object contains a result, i.e. if an algorithm ran and converged.
Should be called before the result is used via `data.v`.
"""
function hasresult(data::PLaplaceData)
    return !ismissing(data.v)
end

"""
$(TYPEDSIGNATURES)

Adds required values from the static data to the output data. In particular
the prolonged boundary and which nodes will not be fixed during the computation.
"""
function add_staticdata!(data::PLaplaceData, S::StaticData)
    data.g = S.g
    data.calculation_nodes = S.calculation_nodes
end

"""
$(TYPEDSIGNATURES)

Adds required values from the iteration data to the output data, especially the result.
"""
function add_algorithmdata!(data::PLaplaceData, algorithm::AlgorithmData)
    data.tsetup = algorithm.tsetup
    data.taux = algorithm.taux
    data.tmain = algorithm.tmain

    data.Naux = algorithm.Naux
    data.Nmain = algorithm.Nmain

    data.msg = algorithm.msg

    data.u = algorithm.solution
    assemble_result!(data)
end

"""
$(TYPEDSIGNATURES)

Writes result vector to a VTK file with the given name.
The file itself will be created, but the path has to exist prior. 
"""
function write_result_to_vtk(filename::String, data::PLaplaceData)
    if !hasresult(data)
        return
    end

    write_to_vtk(data.v, data.mesh, "v", filename, qdim=data.qdim)
end

"""
$(TYPEDSIGNATURES)

Writes result vector to a VTK file with the given name.
The file itself will be created, but the path has to exist prior. 
"""
function write_result_to_txt(filename::String, data::PLaplaceData)
    if !hasresult(data)
        return
    end

    write_to_txt(data.v, data.mesh, filename, qdim=data.qdim)
end

"""
$(TYPEDSIGNATURES)
    
Prints problem data that had to be specified for the algorithm ot the console. 
"""
function print_defaultdata(data::PLaplaceData)
    log_defaultdata(
        LogData(),
        data.mesh,
        data.eps,
        data.p,
        string(data.stepsize)
    )
end

"""
$(TYPEDSIGNATURES)
    
Prints statistics of the algorithm run ot the console. 
"""
function print_statistics(data::PLaplaceData)
    log_statistics(
        LogData(),
        data.tsetup,
        data.taux,
        data.tmain,
        data.Naux,
        data.Nmain,
        data.eps,
        data.msg
    )
end

"""
    write_statistics_header(filename::String; guarded::Bool=false)
    
Clears and writes a header for a statistics log to the given file.
In case the file does not exist, it will be created, but only if the path exists.
If guarded checks before if file already contains a header
and then does not overwrite potentially previous results.
"""
function write_statistics_header(filename::String; guarded::Bool=false)
    fn = occursin(".", filename) ? filename : filename * ".txt"

    if guarded
        check_statistics_header(fn) && return
    end

    open(fn, "w") do file
        write(file, rpad("p",6), "|")
        write(file, rpad("eps",13), "|")
        write(file, rpad("n",7), "|")
        write(file, rpad("m",7), "|")
        write(file, rpad("Scheme",8), "|")
        write(file, rpad("Naux",6), "|")
        write(file, rpad("Nmain",6), "|")
        write(file, rpad("N",7), "|")
        write(file, rpad("t setup",9), "|")
        write(file, rpad("t aux",9), "|")
        write(file, rpad("t main",9), "|")
        write(file, rpad("t sum",9), "|")
        write(file, rpad("Message",10), "\n")
        write(file, repeat("-", 115), "\n")
        write(file, "\$Simulations", "\n")
    end
end

"""
$(TYPEDSIGNATURES)
    
Checks if given file exists and alredy contains a statistics header. 
"""
function check_statistics_header(filename::String) :: Bool
    fn = occursin(".", filename) ? filename : filename * ".txt"

    if !isfile(fn)
        return false
    end

    f = open(fn)
    l = readline(f)
    close(f)
    a = split(l, "|")
    
    !contains(a[1],"p") && return false
    !contains(a[2],"eps") && return false
    !contains(a[13],"Message") && return false

    return true
end

"""
$(TYPEDSIGNATURES)
    
Writes statistics line corresponding to the header to a given log file. 
"""
function write_statistics(filename::String, data::PLaplaceData)
    fn = occursin(".", filename) ? filename : filename * ".txt"

    sp = @sprintf("%06.3f", data.p)
    se = @sprintf("%.7e", data.eps)
    sn = @sprintf("%.07i", data.mesh.nnodes)
    sm = @sprintf("%.07i", data.mesh.nelems)

    st = replace(@sprintf("%-8s", data.stepsize)," "=>"~")
    na = ismissing(data.Naux) ? repeat("~", 6) : @sprintf("%+.05i", data.Naux)
    nm = ismissing(data.Nmain) ? repeat("~", 6) : @sprintf("%+.05i", data.Nmain)
    
    Nc::Int64 = 0
    Nc += ismissing(data.Naux) ? 0 : data.Naux
    if !ismissing(data.Nmain)
        Nc += abs(data.Nmain)
        Nc *= sign(data.Nmain)
    end
    nc = iszero(Nc) ? repeat("~", 7) : @sprintf("%+.06i", Nc)

    ts = ismissing(data.tsetup) ? "~~~~~.~~~" : @sprintf("%09.3f", data.tsetup)
    ta = ismissing(data.taux) ? "~~~~~.~~~" : @sprintf("%09.3f", data.taux)
    tm = ismissing(data.tmain) ? "~~~~~.~~~" : @sprintf("%09.3f", data.tmain)

    tsum::Float64 = 0
    tsum += ismissing(data.tsetup) ? 0 : data.tsetup
    tsum += ismissing(data.taux) ? 0 : data.taux
    tsum += ismissing(data.tmain) ? 0 : data.tmain
    tc = iszero(tsum) ? "~~~~.~~~" : @sprintf("%09.3f", tsum)

    open(fn, "a") do file
        write(file, sp, " ")
        write(file, se, " ")
        write(file, sn, " ")
        write(file, sm, " ")
        write(file, st, " ")
        write(file, na, " ")
        write(file, nm, " ")
        write(file, nc, " ")
        write(file, ts, " ")
        write(file, ta, " ")
        write(file, tm, " ")
        write(file, tc, " ")
        write(file, data.msg, "\n")
    end
end
