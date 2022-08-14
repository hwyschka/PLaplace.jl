using PLaplace, MinFEM, Plots

function Testcase()
    outputPath = "results/SourceProblem1D/"
    mkpath(outputPath)
    write_log_header(outputPath * "log")

    g(x) = 0
    f(x) = 1
    
    mesh = unit_interval(1001)
    dirichletBoundary = select_boundaries(mesh, 1001, 1002)

    ps = [2.0, 3.0, 5.0, 8.0, 15.0, 25.0, 50.0, 100.0]

    for p in ps
        data::PLaplaceData = solve_plaplace(p, mesh, g, dirichletBoundary, f=f)
    
        write_log(outputPath * "log", data)
        write_result_to_vtk(outputPath * "SourceProblem1D_p=$p", data)
        write_to_txt(data.v, mesh, outputPath * "SourceProblem1D_p=$p")
    end

    plt = plot(getindex.(mesh.Nodes,1), evaluate_mesh_function(mesh,f), label="f", c=:green, line=:dash, aspect_ratio=:equal, xlims=(0,1), ylims=(0,0.5))
    
    for p in ps
        N, v = read_from_txt(outputPath * "SourceProblem1D_p=$p.txt")

        plot!(plt, getindex.(N,1), v, label="p=$(trunc(Int, p))")
    end
    
    savefig(plt, outputPath * "1DSouceProblemPlot.svg")
end

Testcase()
