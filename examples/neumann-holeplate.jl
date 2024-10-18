using PLaplace

function testcase()
    output_path = "results/neumann-holeplate/"
    mkpath(output_path)

    statistics_file = output_path * "statistics.txt"
    write_statistics_header(statistics_file, guarded = true)

    mesh = import_mesh("meshes/rectangle_hole.msh")
    dirichlet_boundary = select_boundaries(mesh, 1001, 1002, 1003, 1004)
    neumann_boundary = select_boundaries(mesh, 1005)

    g(x) = 0
    h(x) = x[1]-1

    p::Float64 = 5

    data = solve_plaplace(
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

testcase()
