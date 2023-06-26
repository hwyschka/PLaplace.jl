using PLaplaceDev, MinFEMDev

function testcase()
    outputPath::String = "results/manufactured_solutions/dirichlet_scalar/"
    mkpath(outputPath)

    log_file = outputPath * "log.txt"
    if !isfile(log_file)
        write_log_header(log_file)
    end

    error_file = outputPath * "error.txt"
    if !isfile(error_file)
        write_error_header(error_file)
    end

    eps::Float64 = 1e-6
    stepsize::Stepsize = ADAPTIVE
    kappa::Float64 = 10
    maxIterations::Int64 = 10000
    p::Float64 = 5
    gridpointsPerDim::Int64 = 200

    f(x) = -p * (x[1]^2 + x[2]^2)^((p-2)/2)
    g(x) = 0.5 * (x[1]^2 + x[2]^2)

    mesh::Mesh = unit_square(gridpointsPerDim)
    dirichletBoundary::Set{Boundary} = select_boundaries(mesh)

    outputData = solve_plaplace(
        p,
        mesh,
        g,
        dirichletBoundary,
        f=f, 
        eps=eps,
        stepsize=stepsize,
        kappa=kappa,
        maxIter=maxIterations
    )

    mansol = evaluate_mesh_function(mesh, g)
    
    write_to_vtk(
        [
            outputData.v,
            mansol,
            outputData.v - mansol
        ],
        mesh, 
        ["v", "v_analytic", "v - v_analytic"],
        outputPath * "result_p=$p-n=$gridpointsPerDim-$eps-$stepsize"
    )

    write_log(log_file, outputData)
    write_error(error_file, outputData, compute_errors(mansol, outputData))
end

testcase()
