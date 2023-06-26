"""
    @enum Stepsize

Variable used to specify the stepsize of path following schemes.
"""
@enum Stepsize begin
    SHORT
    LONG
    ADAPTIVE
end 

"""
    mutable struct PLaplaceData

Stores all parameters and results obtained by the solution of a pLaplace problem.

# Fields
- `p::Float64`: problem parameter.
- `qdim::Int64`: image dimension of functions.
- `eps::Float64`: final accuracy of main path-following.
- `stepsize::Stepsize`: stepsize used for the path-following schemes.
- `mesh::Mesh`: FEM mesh of the domain.
- `calculationNodes::Set{Int64}`: nodes of the mesh that have actually been computed and
                                    were not fixed by a Dirichlet boundary.
- `Naux::Int64`: required iterations of the auxilliary path-following
- `Nmain::Int64`: required iterations of the main path-following.
- `taux::Float64`: required time for the auxilliary path-following in seconds.
- `g::Function`: prolonged Dirichlet boundary condition.
- `u::Function`: solution obtained by the algorithm (shifted to zero trace).
- `v::Function`: combined solution to the p-Poisson problem.
- `msg::String`: message if path-following stopped early.

"""
mutable struct PLaplaceData
    "computational domain"
    mesh::Mesh

    "Dirichlet boundary of the computational domain"
    dirichletBoundary::Set{Boundary}

    "Neumann boundary of the computational domain"
    neumannBoundary::Set{Boundary}

    "domain source term"
    f::AbstractVector{Float64}

    "Neumann boundary condition"
    h::AbstractVector{Float64}
    
    "image dimension of solution"
    qdim::Int64

    "PDE parameter"
    p::Float64

    "accuracy"
    eps::Float64

    "stepsize of the interior point method"
    stepsize::Stepsize
    
    "required iterations for the auxilliary path-following"
    Naux::Int64

    "required iterations for the main path-following"
    Nmain::Int64

    "required time for the setup"
    tsetup::Float64

    "required time for the auxilliary path-following"
    taux::Float64

    "required time for the main path-following"
    tmain::Float64

    "nodes on which the solution was computed, i.e. no zero Dirichlet nodes"
    calculationNodes::Set{Int64}

    "Dirichlet boundary condition prolonged to the domain"
    g::AbstractVector{Float64}

    "shifted solution on calculation nodes"
    u::AbstractVector{Float64}

    "solution to the pLaplace problem"
    v::AbstractVector{Float64}

    "notifications from the iteration"
    msg::String
end

PLaplaceData(
    mesh::Mesh,
    dirichletBoundary::Set{Boundary},
    neumannBoundary::Set{Boundary},
    f::AbstractVector{Float64},
    h::AbstractVector{Float64},
    qdim::Int64,
    p::Float64,
    eps::Float64,
    stepsize::Stepsize; 
    Naux::Int64=0,
    Nmain::Int64=0,
    tsetup::Float64=0.0,
    taux::Float64=0.0,
    tmain::Float64=0.0, 
    calculationNodes::Set{Int64}=Set{Int64}(),
    g::AbstractVector{Float64}=Vector{Float64}(), 
    u::AbstractVector{Float64}=Vector{Float64}(),
    v::AbstractVector{Float64}=Vector{Float64}(),
    msg::String=""
) = PLaplaceData(
    mesh,
    dirichletBoundary,
    neumannBoundary,
    f,
    h,
    qdim,
    p,
    eps,
    stepsize,
    Naux,
    Nmain,
    tsetup,
    taux,
    tmain,
    calculationNodes,
    g,
    u,
    v,
    msg
)

mutable struct ConditionData
    "iteration ids in auxilliary path-following"
    AuxilliaryIterations::Array{Int64,1}
                        
    "iteration ids in main path-following"
    MainIterations::Array{Int64,1}
                        
    "condtion value in auxilliary iteration"
    AuxilliaryConditions::Array{Float64,1}
                        
    "condition value in main iteration"
    MainConditions::Array{Float64,1}                    
end

mutable struct ErrorData
    "PDE parameters"
    p::Array{Float64,1}
                        
    "number of gridpoints"
    n::Array{Int64,1}

    "accuracy"
    eps::Array{Float64,1}
                        
    "errors"
    errors::Array{Array{Float64,1},1}
end

"""
    assemble_result!(out::PLaplaceData)

Assembles result vector v = u + g for truncated u vector u_I.
"""
function assemble_result!(out::PLaplaceData)
    n = out.mesh.nnodes
    v = zeros(out.qdim*n)
    k = 1
    for i = 1:n
        v[out.qdim*(i-1)+1:out.qdim*i] = out.g[out.qdim*(i-1)+1:out.qdim*i]
        
        if i in out.calculationNodes
            v[out.qdim*(i-1)+1:out.qdim*i] += out.u[out.qdim*(k-1)+1:out.qdim*k]
            k += 1
        end
    end
    
    out.v = v
end

"""
    add_iterationdata!(data::PLaplaceData, I::IterationData, S::StaticData)

Adds required values from the iteration data to the output data, especially the result.
"""
function add_performance!(outputData::PLaplaceData, 
                            times::Vector{Float64}, iterations::Vector{Int64})
    outputData.tsetup = times[1]
    outputData.taux = times[2]
    outputData.tmain = times[3]

    outputData.Naux = iterations[1]
    outputData.Nmain = iterations[2]
end

"""
    write_result_to_vtk(file_name::String, out::PLaplaceData)

Writes result vector to a VTK file. 
"""
function write_result_to_vtk(file_name::String, out::PLaplaceData)
    write_to_vtk(out.v, out.mesh, "v", file_name, qdim=out.qdim)
end

"""
    write_log_header(file_name::String)
    
Clears and writes a header for a log to the given file. 
"""
function write_log_header(file_name::String)
    fn = occursin(".", file_name) ? file_name : file_name * ".txt"

    open(fn, "w") do file
        write(file, rpad("p",6), "|")
        write(file, rpad("eps",13), "|")
        write(file, rpad("n",7), "|")
        write(file, rpad("m",7), "|")
        write(file, rpad("Scheme",8), "|")
        write(file, rpad("Naux",5), "|")
        write(file, rpad("Nmain",5), "|")
        write(file, rpad("N",6), "|")
        write(file, rpad("t setup",9), "|")
        write(file, rpad("t aux",9), "|")
        write(file, rpad("t main",9), "|")
        write(file, rpad("t sum",9), "|")
        write(file, rpad("Message",10), "\n")
        write(file,repeat("-", 115), "\n")
    end
end

"""
    write_log(file_name::String, outputData::PLaplaceData)
    
Writes data corresponding to the header to a given log file. 
"""
function write_log(file_name::String, outputData::PLaplaceData)
    fn = occursin(".", file_name) ? file_name : file_name * ".txt"

    sp = @sprintf("%06.3f", outputData.p)
    se = @sprintf("%.7e", outputData.eps)
    sn = @sprintf("%.07i", outputData.mesh.nnodes)
    sm = @sprintf("%.07i", outputData.mesh.nelems)

    st = replace(@sprintf("%-8s", outputData.stepsize)," "=>"~")
    na = @sprintf("%.05i", outputData.Naux)
    nm = @sprintf("%.05i", outputData.Nmain)
    nc = @sprintf("%.06i", outputData.Naux + outputData.Nmain)

    ts = @sprintf("%09.3f", outputData.tsetup)
    ta = @sprintf("%09.3f", outputData.taux)
    tm = @sprintf("%09.3f", outputData.tmain)
    sum = outputData.tsetup + outputData.taux + outputData.tmain
    tc = @sprintf("%09.3f", sum)

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
        write(file, outputData.msg, "\n")
    end
end

"""
    write_error_header(file_name::String)
    
Clears and writes a header for an error log to the given file. 
"""
function write_error_header(file_name::String)
    fn = occursin(".", file_name) ? file_name : file_name * ".txt"

    open(fn, "w") do file
        write(file, rpad("p",6), "|")
        write(file, rpad("n",7), "|")
        write(file, rpad("eps",13), "|")
        write(file, rpad("J-error",14), "|")
        write(file, rpad("L1-error",13), "|")
        write(file, rpad("L2-error",13), "|")
        write(file, rpad("LInf-error",13), "\n")
        write(file, repeat("-", 85), "\n")
        write(file, "\$Simulations", "\n")
    end
end

"""
    write_error(file_name::String, outputData::PLaplaceData, error::Float64)
    
Writes a given error along other information corresponding to the header to a given file.
"""
function write_error(file_name::String, outputData::PLaplaceData, error::Array{Float64,1})
    fn = occursin(".", file_name) ? file_name : file_name * ".txt"
    
    sp = @sprintf("%06.3f", outputData.p)
    sn = @sprintf("%.07i", outputData.mesh.nnodes)
    ee = @sprintf("%.7e", outputData.eps)
    oe = @sprintf("%+.7e", error[1])
    e1 = @sprintf("%.7e", error[2])
    e2 = @sprintf("%.7e", error[3])
    ei = @sprintf("%.7e", error[4])
    open(fn, "a") do file
        write(file, sp, " ")
        write(file, sn, " ")
        write(file, ee, " ")
        write(file, oe, " ")
        write(file, e1, " ")
        write(file, e2, " ")
        write(file, ei, "\n")
    end
end

"""
    write_step_header(fileName::String, scheme::String)
    
Clears and writes a header for an stepwise log to the given file. 
"""
function write_step_header(fileName::String, scheme::String)
    open(fileName * ".txt", "w") do file
        write(file,"Scheme: ", scheme, "\n")
        write(file, rpad("Step",7), "|")
        write(file, rpad("Value",9), "|")
        write(file, rpad("n",7), "|")
        write(file, rpad("L2-error",12), "\n")
        write(file,repeat("-", 38), "\n")
    end
end

"""
    write_step(fileName::String, iter::Int64, value::Float64)
    
Writes a given value corresponding to an iteration count to a given file. 
"""
function write_step(fileName::String, iter::Int64, value::Float64)
    id = @sprintf("%.04i", iter)
    val = @sprintf("%.7e", value)
    open(fileName * ".txt", "a") do file
        write(file, id, " ")
        write(file, string(val), "\n")
    end
end

"""
    read_step(fileName::String)
    
Returns iterations and values from a file written by write_step. 
"""
function read_stepfile(fileName::String)
    totalIterations = 0
    open(fileName) do f
        while (!eof(f))
            l = readline(f)
            totalIterations += 1
        end
    end

    iteration = Array{Int64,1}(undef,totalIterations)
    value = Array{Float64,1}(undef,totalIterations)
    
    f = open(fileName)

    k = 1
    while (!eof(f))
        l = readline(f)
        a = split(l, " ")
        
        iteration[k] = parse(Int64, a[1])
        value[k] = parse(Float64, a[2])
        
        k += 1
    end

    close(f)

    return iteration, value
end

"""
    read_condition(fileName::String)
    
Returns condition data by reading the two condition files generated by solve_plaplace
when handed over the given fileName (i.e. without _aux.txt or _main.txt ending). 
"""
function read_condition(fileName::String)
    iterations_aux, conditions_aux = read_stepfile(fileName * "_aux.txt")
    iterations_main, conditions_main = read_stepfile(fileName * "_main.txt")
    
    
    return ConditionData(iterations_aux, iterations_main,
                            conditions_aux, conditions_main)
end

"""
    read_error(fileName::String)
    
Returns error data when handed an file generated by PLaplace. 
"""
function read_error(fileName::String)
    p = Array{Float64,1}()
    n = Array{Int64,1}()
    eps = Array{Float64,1}()
    errors = Array{Array{Float64,1},1}()
    
    f = open(fileName)

    while (!eof(f) && (l = readline(f)) != "\$Simulations")
    end

    while (!eof(f))
        l = readline(f)
        a = split(l, " ")
        
        push!(p, parse(Float64, a[1]))
        push!(n, parse(Int64, a[2]))
        push!(eps, parse(Float64, a[3]))
        
        e = [
            parse(Float64, a[4]),
            parse(Float64, a[5]),
            parse(Float64, a[6]),
            parse(Float64, a[7])
        ]
        push!(errors, e)      
    end

    close(f)
    
    return ErrorData(p, n, eps, errors)
end

function append!(A::ErrorData, B::ErrorData)
    append!(A.p, B.p)
    append!(A.n, B.n)
    append!(A.eps, B.eps)
    append!(A.errors, B.errors)
end

"""
    print_defaultdata(outputData::PLaplaceData)
    
Prints problem data that had to be specified for the algorithm ot the console. 
"""
function print_defaultdata(outputData::PLaplaceData)
    println("========================================================")
    println("Mesh: #nodes = ", outputData.mesh.nnodes)
    println("      #elements = ", outputData.mesh.nelems)
    println("Accuracy: eps = ", outputData.eps)
    println("PDE Parameter: p = ", outputData.p)
    println("Stepsize: ", outputData.stepsize)
    println("========================================================")
end

"""
    print_statistics(outputData::PLaplaceData)
    
Prints required time and iterations of the solution process to the console.
"""
function print_statistics(outputData::PLaplaceData)
    tsetup = @sprintf("%.3f", outputData.tsetup)
    taux = @sprintf("%.3f", outputData.taux)
    tmain = @sprintf("%.3f", outputData.tmain)
    tsum = @sprintf("%.3f", outputData.taux + outputData.tmain + outputData.tsetup)
    N = outputData.Naux + outputData.Nmain
    println("========================================================")
    println("Calculation statistics")
    println("========================================================")
    println("Setup in $tsetup seconds.")
    println("$(outputData.Naux) auxilliary iterations in $taux seconds.")
    println("$(outputData.Nmain) main iterations in $tmain seconds.")
    println("Setup and $N iterations combined in $tsum seconds.")
    println("Stop Message: $(outputData.msg)")
end
