using PLaplace, MinFEM

function Testcase()
    outputPath = "results/NeumannProblem2DVector/"
    mkpath(outputPath)

    g(x) = [0,0]
    h(x) = (x[1]^2 + x[2]^2) * [1,1]

    mesh = import_mesh("meshes/square.msh")
    dirichletBoundary = select_boundaries(mesh, 1001, 1004)
    neumannBoundary = select_boundaries(mesh, 1002, 1003)

    p::Float64 = 5

    data::PLaplaceData = solve_plaplace(p, mesh, g, dirichletBoundary, neumannBoundary=neumannBoundary, h=h, qdim=2)
    
    write_result_to_vtk(outputPath * "NeumannProblem2DVector_p=$p", data)
end

Testcase()