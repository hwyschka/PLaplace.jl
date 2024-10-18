using PLaplace

function testcase()
    output_path = "results/dirichlet-square/"
    mkpath(output_path)
    
    statistics_file = output_path * "statistics.txt"
    write_statistics_header(statistics_file, guarded = true)

    mesh = import_mesh("meshes/square.msh")
    dirichlet_boundary = select_boundaries(mesh)

    g(x) = x[1]^2

    for p in [1.0, 2.0, 3.0, 5.0, 8.0, 15.0]
        data = solve_plaplace(p, mesh, g, dirichlet_boundary)
    
        write_statistics(statistics_file, data)
        write_result_to_vtk(output_path * "result_p=$p", data)
    end
end

testcase()
