"""
    PLaplace

A Solver for Problems of the p-Laplacian PDE Operator.

This package imports the following packages:
$(IMPORTS)
"""
module PLaplace
    using MinFEM
    using LinearAlgebra, SparseArrays, InvertedIndices
    using IncompleteLU, IterativeSolvers, Preconditioners
    using Printf, WriteVTK
    using DocStringExtensions

    include("utility.jl")

    include("linearsystems.jl")
    include("staticdata.jl")


    include("barrierfunction.jl")
    include("barriers/finite.jl")

    include("logbase.jl")
    include("logdata.jl")
    include("algorithmlog.jl")

    include("tracking.jl")
    include("iterationdata.jl")
    include("algorithmdata.jl")
    include("pathfollowing.jl")

    include("plaplacedata.jl")

    include("problem.jl")
    include("errordata.jl")

    include("algorithm.jl")

    export  PLaplaceData,
            AlgorithmLogData,
            ErrorData

    export  objective_functional,
            compute_errors,
            write_error_header,
            check_error_header,
            write_error,
            read_error,
            append!

    export  solve_plaplace

    export  Mesh,
            Boundary,
            import_mesh,
            select_boundaries

    export  hasresult,
            write_result_to_vtk,
            write_result_to_txt,
            print_defaultdata,
            print_statistics,
            write_statistics_header,
            check_statistics_header,
            write_statistics,
            read_algorithmlog

    export  assemble_derivativetensor,
            assemble_derivativetensor_boundary,
            assemble_derviativetensor_modified,
            compute_derivative,
            compute_normalderivative,
            xpnorm
        
    export  compute_prolongation_harmonic,
            compute_prolongation_zero

    export  compute_lipschitzconstant_boundary

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
