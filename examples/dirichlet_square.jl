using PLaplace, MinFEM

function Testcase()
    outputPath = "results/DirichletProblem2D/"
    mkpath(outputPath)
    write_log_header(outputPath * "log")

    g(x) = x[1]^2

    mesh = import_mesh("meshes/square.msh")
    dirichletBoundary = select_boundaries(mesh)

    for p in [2.0, 3.0, 5.0, 8.0, 15.0, 25.0]
        data::PLaplaceData = solve_plaplace(p, mesh, g, dirichletBoundary)
    
        write_result_to_vtk(outputPath * "DirichletProblem2D_p=$p", data)
        write_log(outputPath * "log", data)
    end
end

Testcase()