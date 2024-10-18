using PLaplace, MinFEM

function testcase()
    # Setup output
    output_path::String = "results/manufactured_solutions/neumann_scalar/"
    mkpath(output_path)

    statistics_file = output_path * "statistics.txt"
    write_statistics_header(statistics_file, guarded = true)

    error_file = output_path * "errors.txt"
    write_error_header(error_file, guarded=true)

    # Select PDE parameter(s)
    #p::Float64 = 2.0
    for p::Float64 in [2, 2.1, 2.5, 3, 4, 4.1, 4.8, 5, 6, 7, 8, 9, 10, 11, 12]

    f(x) = -p * (x[1]^2 + x[2]^2)^((p-2)/2)
    h(x) = x[1] > 0 ? 
        x[2] * (x[1]^2 + x[2]^2)^((p-2)/2) :
        -x[1] * (x[1]^2 + x[2]^2)^((p-2)/2)
    g(x) = 0.5 * (x[1]^2 + x[2]^2)

    # Select number(s) of gridpoints [ignored for gmsh]
    ndim::Int64 = 200
    #for ndim in [10, 20, 50, 100, 200, 500]

    # Select unit square mesh or custom mesh
    mesh::Mesh = unit_square(ndim)
    #mesh::Mesh = import_mesh("meshes/square.msh")

    # Define boundary sets
    dirichlet_boundary::Set{Boundary} = select_boundaries(mesh, 1001, 1004)
    neumann_boundary::Set{Boundary} = select_boundaries(mesh, 1002, 1003)

    _h = evaluate_quadrature_function_boundary(mesh, h, neumann_boundary)

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
        h = _h,
        neumann_boundary = neumann_boundary,
        eps = eps,
        stepsize = stepsize,
        kappa = kappa,
        maxiterations = maxiterations,
        logfile = output_path * "log_p=$p-n=$(mesh.nnodes)-$eps-$stepsize",
        logcondition = true
    )

    write_statistics(statistics_file, data)

    if hasresult(data)
        mansol = evaluate_mesh_function(mesh, g)
        write_error(error_file, data, compute_errors(mansol, data))
    
        write_to_vtk(
            [
                data.v,
                mansol,
                data.v - mansol
            ],
            mesh, 
            ["v", "v_analytic", "v - v_analytic"],
            output_path * "result_p=$p-n=$(mesh.nnodes)-$eps-$stepsize"
        )
    end

    end
    #end
    #end
    #end
end

testcase()
