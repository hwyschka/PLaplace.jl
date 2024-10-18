using PLaplace, MinFEM

function testcase()
    # Setup output
    output_path::String = "results/manufactured_solutions/dirichlet_vector/"
    mkpath(output_path)

    statistics_file = output_path * "statistics.txt"
    write_statistics_header(statistics_file, guarded = true)

    error_file = output_path * "errors.txt"
    write_error_header(error_file, guarded=true)

    # Select PDE parameter(s)
    #p::Float64 = 2.0
    for p::Float64 in [2, 2.1, 2.5, 3, 4, 4.1, 4.8, 5, 6, 7, 8, 9]

    f(x) = -p * 2^((p-2)/2) * (x[1]^2 + x[2]^2)^((p-2)/2) * [1,1]
    g(x) = 0.5 * (x[1]^2 + x[2]^2) * [1,1]

    # Finite manufactured solution II (switched signs in 2nd component)
    #f(x) = -1*2^((3/2)*p - 2) * [p*(x[1]^2 + x[2]^2)^((p-2)/2), (x[1]^2 + x[2]^2)^((p-2)/2 -1)*((p-2)*x[1]^2-p*x[2]^2)]
    #g(x) = [x[1]^2 + x[2]^2, x[1]^2 - x[2]^2]

    # Select number(s) of gridpoints [ignored for gmsh]
    ndim::Int64 = 200
    #for ndim in [10, 20, 50, 100, 200, 500]

    # Select unit square mesh or custom mesh
    mesh::Mesh = unit_square(ndim)
    #mesh::Mesh = import_mesh("meshes/square.msh")

    # Define boundary set
    dirichlet_boundary::Set{Boundary} = select_boundaries(mesh)

    # Select accuracy(s)
    eps::Float64 = 1e-6
    #for eps::Float64 in [1e-5, 1e-6, 1e-7]

    # Select stepsize(s)
    stepsize::Stepsize = ADAPTIVE
    #for stepsize::Stepsize in [SHORT, LONG, ADAPTIVE]
    
    # Adjust other parameters
    kappa::Float64 = 10
    maxiterations::Int64 = 10000


    # Solve p-Laplace equation
    data = solve_plaplace(
        p,
        mesh,
        g,
        dirichlet_boundary,
        f = f,
        qdim = 2, 
        eps = eps,
        stepsize = stepsize,
        kappa = kappa,
        maxiterations = maxiterations,
        logfile = output_path * "log_p=$p-n=$(mesh.nnodes)-$eps-$stepsize",
        logcondition = true
    )

    write_statistics(statistics_file, data)

    if hasresult(data)
        mansol = evaluate_mesh_function(mesh, g, qdim=2)
        write_error(error_file, data, compute_errors(mansol, data))
    
        write_to_vtk(
            [
                data.v,
                mansol,
                data.v - mansol
            ],
            mesh, 
            ["v", "v_analytic", "v - v_analytic"],
            output_path * "result_p=$p-n=$(mesh.nnodes)-$eps-$stepsize",
            qdim = data.qdim
        )
    end

    end
    #end
    #end
    #end
end

testcase()
