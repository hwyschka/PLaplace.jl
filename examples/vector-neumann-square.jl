using PLaplace, MinFEM

function testcase()
    output_path = "results/vector-neumann-square/"
    mkpath(output_path)

    statistics_file = output_path * "statistics.txt"
    write_statistics_header(statistics_file, guarded = true)

    mesh = unit_square(50)
    dirichlet_boundary = select_boundaries(mesh, 1001, 1004)
    neumann_boundary = select_boundaries(mesh, 1002, 1003)

    g(x) = [0,0]
    h(x) = (x[1]^2 + x[2]^2) * [1,1]

    p::Float64 = 5

    data = solve_plaplace(
        p,
        mesh,
        g,
        dirichlet_boundary,
        neumann_boundary = neumann_boundary,
        h = h,
        qdim = 2
    )
    
    write_statistics(statistics_file, data)
    write_result_to_vtk(output_path * "result_p=$p", data)
end

testcase()
