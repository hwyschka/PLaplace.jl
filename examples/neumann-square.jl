using PLaplace

function testcase()
    output_path = "results/neumann-square/"
    mkpath(output_path)
    
    statistics_file = output_path * "statistics.txt"
    write_statistics_header(statistics_file, guarded = true)

    mesh = import_mesh("meshes/square.msh")
    dirichlet_boundary = select_boundaries(mesh, 1001, 1004)
    neumann_boundary = select_boundaries(mesh, 1002, 1003)

    h(x) = x[2]^2 - x[1]^2
    g(x) = 0

    for p in [2.0, 3.0, 5.0, 8.0, 15.0, 25.0]
        data  = solve_plaplace(
            p,
            mesh,
            g,
            dirichlet_boundary,
            h = h,
            neumann_boundary = neumann_boundary
        )
    
        write_statistics(statistics_file, data)
        write_result_to_vtk(output_path * "result_p=$p", data)
    end
end

testcase()
