using PLaplace

function testcase()
    output_path = "results/vector-neumann-cube/"
    mkpath(output_path)

    statistics_file = output_path * "statistics.txt"
    write_statistics_header(statistics_file, guarded = true)

    mesh = import_mesh("meshes/cube.msh")
    neumann_boundary = select_boundaries(mesh, 1003, 1004, 1005)
    dirichlet_boundary = select_boundaries(mesh, 1001, 1002, 1006)

    g(x) = [0,0,0]
    function h(x)
        a = sqrt(sum(x.^2))
        if x[1] > 0
            if x[2] > 0
                return a * [0,0,-1]
            else
                return a * [0,-1,0]
            end
        else
            return a * [-1,0,0]
        end
    end

    p::Float64 = 5

    data  = solve_plaplace(
        p,
        mesh,
        g,
        dirichlet_boundary,
        h = h,
        neumann_boundary = neumann_boundary,
        qdim = 3
    )
    
    write_statistics(statistics_file, data)
    write_result_to_vtk(output_path * "result_p=$p", data)
end

testcase()
