using PLaplace, MinFEM
using LinearAlgebra

function test_scalar(p::Float64)    
    g(x) = 0.5 * norm(x, 2)^2

    f1(x) = -(p-1) * norm(x,2)^(p-2)
    f2(x) = -p * norm(x,2)^(p-2)
    f3(x) = -(p+1) * norm(x,2)^(p-2)

    mesh = import_mesh("test_line.msh")
    dirichletBoundary = select_boundaries(mesh)
    
    outputData = solve_plaplace(p, mesh, g, dirichletBoundary, f=f1, consoleOutput=false)
    pnorm(2.0, outputData.v - evaluate_mesh_function(mesh, g), mesh) > 1e-3 && return false

    mesh = import_mesh("test_square.msh")
    dirichletBoundary = select_boundaries(mesh)
    
    outputData = solve_plaplace(p, mesh, g, dirichletBoundary, f=f2, consoleOutput=false)
    pnorm(2.0, outputData.v - evaluate_mesh_function(mesh, g), mesh) > 1e-3 && return false


    mesh = import_mesh("test_cube.msh")
    dirichletBoundary = select_boundaries(mesh)
    
    outputData = solve_plaplace(p, mesh, g, dirichletBoundary, f=f3, consoleOutput=false)
    pnorm(2.0, outputData.v - evaluate_mesh_function(mesh, g), mesh) > 1e-3 && return false

    return true
end

@test test_scalar(3.0)
@test test_scalar(5.0)
