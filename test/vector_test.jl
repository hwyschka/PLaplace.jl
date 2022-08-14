using PLaplace, MinFEM
using LinearAlgebra

function test_vector(p::Float64)    
    g1(x) = 0.5 * norm(x, 2)^2
    g2(x) = 0.5 * norm(x, 2)^2 * [1,1]
    g3(x) = 0.5 * norm(x, 2)^2 * [1,1,1]

    f1(x) = -(p-1) * norm(x,2)^(p-2)
    f2(x) = -p * 2^((p-2)/2) * norm(x,2)^(p-2) * [1,1]
    f3(x) = -(p+1) * 3^((p-2)/2) * norm(x,2)^(p-2) * [1,1,1]

    f(x) = -p * 2^((p-2)/2) * (x[1]^2 + x[2]^2)^((p-2)/2) * [1,1]

    mesh = import_mesh("test_line.msh")
    dirichletBoundary = select_boundaries(mesh)
    
    outputData = solve_plaplace(p, mesh, g1, dirichletBoundary, f=f1, qdim=1, consoleOutput=false)
    pnorm(2.0, outputData.v - evaluate_mesh_function(mesh, g1, qdim=1), mesh, qdim=1) > 1e-3 && return false

    mesh = import_mesh("test_square.msh")
    dirichletBoundary = select_boundaries(mesh)
    
    outputData = solve_plaplace(p, mesh, g2, dirichletBoundary, f=f2, qdim=2, consoleOutput=false)
    pnorm(2.0, outputData.v - evaluate_mesh_function(mesh, g2, qdim=2), mesh, qdim=2) > 1e-3 && return false


    mesh = import_mesh("test_cube.msh")
    dirichletBoundary = select_boundaries(mesh)
    
    outputData = solve_plaplace(p, mesh, g3, dirichletBoundary, f=f3, qdim=3, consoleOutput=false)
    pnorm(2.0, outputData.v - evaluate_mesh_function(mesh, g3, qdim=3), mesh, qdim=3) > 1e-3 && return false

    return true
end

@test test_vector(3.0)
@test test_vector(5.0)
