using PLaplace, MinFEM

function Testcase()
    outputPath = "results/HolePlateProblem2D/"
    mkpath(outputPath)

    g(x) = 0
    h(x) = x[1]-1

    mesh = import_mesh("meshes/rectangle_hole.msh")
    dirichletBoundary = select_boundaries(mesh, 1001, 1002, 1003, 1004)
    neumannBoundary = select_boundaries(mesh, 1005)

    p::Float64 = 5

    data::PLaplaceData = solve_plaplace(p, mesh, g, dirichletBoundary, h=h, neumannBoundary=neumannBoundary)
    
    write_result_to_vtk(outputPath * "HolePlateProblem2D_p=$p", data)
end

Testcase()
