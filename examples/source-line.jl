using PLaplace, MinFEM, Plots

function testcase()
    output_path = "results/source-line/"
    mkpath(output_path)
    
    statistics_file = output_path * "statistics.txt"
    write_statistics_header(statistics_file, guarded = true)

    mesh = unit_interval(1001)
    dirichlet_boundary = select_boundaries(mesh, 1001, 1002)

    g(x) = 0
    f(x) = 1

    ps = [2.0, 3.0, 5.0, 8.0, 15.0, 25.0, 50.0, 100.0]

    for p in ps
        data::PLaplaceData = solve_plaplace(p, mesh, g, dirichlet_boundary, f=f)
    
        write_log(output_path * "log", data)
        write_result_to_vtk(output_path * "result_p=$p", data)
        write_to_txt(data.v, mesh, output_path * "result_p=$p")
    end

    plt = plot(
        getindex.(mesh.Nodes, 1),
        evaluate_mesh_function(mesh, f),
        label = "f",
        c = :green,
        line = :dash,
        aspect_ratio = :equal,
        xlims = (0,1),
        ylims = (0,0.5)
    )
    
    for p in ps
        N, v = read_from_txt(output_path * "result_p=$p.txt")

        plot!(plt, getindex.(N,1), v, label="p=$(trunc(Int, p))")
    end
    
    savefig(plt, output_path * "plot.svg")
end

testcase()
