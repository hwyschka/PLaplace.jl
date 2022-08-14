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
- `calculationNodes::Set{Int64}`: nodes of the mesh that have actually been computed and were not fixed by a Dirichlet boundary.
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

PLaplaceData(mesh::Mesh, qdim::Int64, p::Float64, eps::Float64, stepsize::Stepsize; 
                Naux::Int64=0, Nmain::Int64=0, tsetup::Float64=0.0, taux::Float64=0.0, tmain::Float64=0.0, 
                calculationNodes::Set{Int64}=Set{Int64}(), g::AbstractVector{Float64}=Vector{Float64}(), 
                u::AbstractVector{Float64}=Vector{Float64}(), v::AbstractVector{Float64}=Vector{Float64}(),
                msg::String="") =
                    PLaplaceData(mesh, qdim, p, eps, stepsize, Naux, Nmain, tsetup, taux, tmain, calculationNodes, g, u, v, msg)

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
function add_performance!(outputData::PLaplaceData, times::Vector{Float64}, iterations::Vector{Int64})
    outputData.tsetup = times[1]
    outputData.taux = times[2]
    outputData.tmain = times[3]

    outputData.Naux = iterations[1]
    outputData.Nmain = iterations[2]
end

"""
    write_result_to_vtk(fileName::String, out::PLaplaceData)

Writes result vector to a VTK file. 
"""
function write_result_to_vtk(fileName::String, out::PLaplaceData)
    write_to_vtk(out.v, out.mesh, "v", fileName, qdim=out.qdim)
end

"""
    write_log_header(fileName::String)
    
Clears and writes a header for a log to the given file. 
"""
function write_log_header(fileName::String)
    open(fileName * ".txt", "w") do file
        write(file,"p      |eps      |n      |m       |Stepsize |Naux  |Nmain |N      |taux     |tmain    |tsetup   | Message\n")
        write(file,"---------------------------------------------------------------------------------------------------------\n")
    end
end

"""
    write_log(fileName::String, outputData::PLaplaceData)
    
Writes data corresponding to the header to a given log file. 
"""
function write_log(fileName::String, outputData::PLaplaceData)
    eps = @sprintf("%.1e", outputData.eps)
    open(fileName * ".txt", "a") do file
        write(file, rpad(string(outputData.p),7), " ", lpad(string(eps),9), " ", lpad(string(outputData.mesh.nnodes),7), " ", lpad(string(outputData.mesh.nelems),8), " ", lpad(string(outputData.stepsize),9), " ", lpad(string(outputData.Naux),6), " ", lpad(string(outputData.Nmain),6), " ", lpad(string(outputData.Nmain+outputData.Naux),7), " ", lpad(string(round(outputData.taux, digits=3)),9), " ", lpad(string(round(outputData.tmain, digits=3)),9), " ", lpad(string(round(outputData.tsetup, digits=3)),9), " ", outputData.msg, "\n")
    end
end

"""
    write_error_header(fileName::String)
    
Clears and writes a header for an error log to the given file. 
"""
function write_error_header(fileName::String)
    open(fileName * ".txt", "w") do file
        write(file,"p      |eps      |n      |L2-error    \n")
        write(file,"--------------------------------------\n")
    end
end

"""
    write_error(fileName::String, outputData::PLaplaceData, error::Float64)
    
Writes a given error along other information corresponding to the header to a given file. 
"""
function write_error(fileName::String, outputData::PLaplaceData, error::Float64)
    err = @sprintf("%.5e", error)
    open(fileName * ".txt", "a") do file
        write(file, lpad(string(outputData.p),7), " ", lpad(string(outputData.eps),9), " ", lpad(string(outputData.mesh.nnodes),7), " ", lpad(string(err),12), "\n")
    end
end

"""
    write_step_header(fileName::String, scheme::String)
    
Clears and writes a header for an stepwise log to the given file. 
"""
function write_step_header(fileName::String, scheme::String)
    open(fileName * ".txt", "w") do file
        write(file,"Scheme: ", scheme, "\n")
        write(file,"step  |val      |n      |L2-error    \n")
        write(file,"----------------------------------\n")
    end
end

"""
    write_step(fileName::String, iter::Int64, value::Float64)
    
Writes a given value corresponding to an iteration count to a given file. 
"""
function write_step(fileName::String, iter::Int64, value::Float64)
    val = @sprintf("%.7e", value)
    open(fileName * ".txt", "a") do file
        write(file, string(iter), ", ", string(val), "\n")
    end
end

"""
    print_defaultdata(outputData::PLaplaceData)
    
Prints problem data that had to be specified for the algorithm ot the console. 
"""
function print_defaultdata(outputData::PLaplaceData)
    println("========================================================")
    println("Mesh: #nodes = ", outputData.mesh.nnodes,", #elements = ", outputData.mesh.nelems)
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
    println("========================================================")
    println("Calculation statistics")
    println("========================================================")
    println("Setup in $tsetup seconds.")
    println("$(outputData.Naux) auxilliary iterations in $taux seconds.")
    println("$(outputData.Nmain) main iterations in $tmain seconds.")
    println("Setup and $(outputData.Naux + outputData.Nmain) iterations combined in $tsum seconds.")
    println("Stop Message: $(outputData.msg)")
end
