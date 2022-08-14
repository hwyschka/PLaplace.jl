module PLaplaceDev

using MinFEM
using LinearAlgebra, SparseArrays, InvertedIndices
using IncompleteLU, IterativeSolvers, Preconditioners
using Printf, WriteVTK

include("utility.jl")
include("output.jl")
include("linearsystems.jl")
include("staticdata.jl")
include("iterationdata.jl")
include("pathfollowing.jl")
include("algorithm.jl")

export  PLaplaceData

export  solve_plaplace

export  Mesh,
        Boundary,
        import_mesh,
        select_boundaries

export  write_result_to_vtk,
        write_log_header,
        write_log,
        write_error_header,
        write_error,
        print_defaultdata,
        print_statistics

export  assemble_derivativetensor,
        assemble_derivativetensor_boundary,
        assemble_derviativetensor_modified,
        compute_derivative,
        compute_normalderivative,
        xpnorm
        
export  compute_prolongation_harmonic,
        compute_prolongation_zero

export Stepsize
for s in instances(Stepsize)
    @eval export $(Symbol(s))
end 

export LinearSolver
for s in instances(LinearSolver)
    @eval export $(Symbol(s))
end 

export Preconditioner
for s in instances(Preconditioner)
    @eval export $(Symbol(s))
end

end
