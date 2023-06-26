using PLaplace, MinFEM

function Testcase()
    outputPath = "results/NeumannProblem2D/"
    mkpath(outputPath)
    write_log_header(outputPath * "log")

    h(x) = x[1]^1 - x[2]^2
    g(x) = 0

    mesh = import_mesh("../meshes/square.msh")
    dirichletBoundary = select_boundaries(mesh, 1001, 1004)
    neumannBoundary = select_boundaries(mesh, 1002, 1003)

    for p in [2.0, 3.0, 5.0, 8.0, 15.0, 25.0]
        data::PLaplaceData  = solve_plaplace(
            p,
            mesh,
            g,
            dirichletBoundary,
            h=h,
            neumannBoundary=neumannBoundary
        )
    
        write_result_to_vtk(outputPath * "NeumannProblem2D_p=$p", data)
        write_log(outputPath * "log", data)
    end
end

Testcase()
