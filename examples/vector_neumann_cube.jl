using PLaplace, MinFEM

function Testcase()
    outputPath = "results/NeumannProblem3DVector/"
    mkpath(outputPath)

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

    mesh = import_mesh("../meshes/cube.msh")
    neumannBoundary = select_boundaries(mesh, 1003, 1004, 1005)
    dirichletBoundary = select_boundaries(mesh, 1001, 1002, 1006)

    p::Float64 = 5

    data::PLaplaceData  = solve_plaplace(
        p,
        mesh,
        g,
        dirichletBoundary,
        h=h,
        neumannBoundary=neumannBoundary,
        qdim=3
    )
    
    write_result_to_vtk(outputPath * "NeumannProblem3DVector_p=$p", data)
end

Testcase()
