using PLaplace, MinFEM

function test_utility()
    f(x) = -x[1]
    g(x) = [x[1]+x[2]+x[3]; x[1]-x[2]+x[3]; x[1]+x[2]-x[3]]
    h(x) = 0.5 .* [x[1], -x[3], x[2]]

    mesh::Mesh = import_mesh("test_line.msh")
    fh = evaluate_mesh_function(mesh, f)

    d1 = -ones(Float64, mesh.nelems)

    D = assemble_derivativetensor(mesh)
    Df = compute_derivative(D, fh)
    any(abs.(Df[(1,1)] .- d1) .> 1e-15) && return false

    Dnf = compute_normalderivative(mesh, select_boundaries(mesh), 0.5 .* fh)
    any(abs.(Dnf[1] .- [0.5, -0.5]) .> 1e-15) && return false


    mesh = import_mesh("test_square.msh")
    fh = evaluate_mesh_function(mesh, f)
    
    d1 = -ones(Float64, mesh.nelems)
    d2 = zeros(Float64, mesh.nelems)

    D = assemble_derivativetensor(mesh)
    Df = compute_derivative(D, fh)
    any(abs.(Df[(1,1)] .- d1) .> 1e-14) && return false
    any(abs.(Df[(2,1)] .- d2) .> 1e-14) && return false

    solDnf = zeros(mesh.nboundelems)
    for bel in extract_elements(select_boundaries(mesh, 1003))
        solDnf[bel] = 0.5
    end
    for bel in extract_elements(select_boundaries(mesh, 1004))
        solDnf[bel] = -0.5
    end

    Dnf = compute_normalderivative(mesh, select_boundaries(mesh), 0.5 .* fh)
    any(abs.(Dnf[1] .- solDnf) .> 1e-14) && return false


    mesh = import_mesh("test_cube.msh")
    fh = evaluate_mesh_function(mesh, f)
    gh = evaluate_mesh_function(mesh, g, qdim=3)
    hh = evaluate_mesh_function(mesh, h, qdim=3)

    d1 = -ones(Float64, mesh.nelems)
    d2 = zeros(Float64, mesh.nelems)

    D = assemble_derivativetensor(mesh)
    Df = compute_derivative(D, fh)
    any(abs.(Df[(1,1)] .- d1) .> 1e-14) && return false
    any(abs.(Df[(2,1)] .- d2) .> 1e-14) && return false
    any(abs.(Df[(3,1)] .- d2) .> 1e-14) && return false 
    
    d1 = -ones(Float64, mesh.nelems)
    d2 = ones(Float64, mesh.nelems)

    D = assemble_derivativetensor(mesh, qdim=3)
    Dg = compute_derivative(D, gh)
    any(abs.(Dg[(1,1)] .- d2) .> 1e-14) && return false
    any(abs.(Dg[(1,2)] .- d2) .> 1e-14) && return false
    any(abs.(Dg[(1,3)] .- d2) .> 1e-14) && return false
    any(abs.(Dg[(2,1)] .- d2) .> 1e-14) && return false
    any(abs.(Dg[(2,2)] .- d1) .> 1e-14) && return false
    any(abs.(Dg[(2,3)] .- d2) .> 1e-14) && return false
    any(abs.(Dg[(3,1)] .- d2) .> 1e-14) && return false
    any(abs.(Dg[(3,2)] .- d2) .> 1e-14) && return false
    any(abs.(Dg[(3,3)] .- d1) .> 1e-14) && return false

    solDnf = zeros(mesh.nboundelems)
    solDnh1 = zeros(mesh.nboundelems)
    solDnh2 = zeros(mesh.nboundelems)
    solDnh3 = zeros(mesh.nboundelems)
    for bel in extract_elements(select_boundaries(mesh, 1001))
        solDnh2[bel] = -0.5
    end
    for bel in extract_elements(select_boundaries(mesh, 1002))
        solDnf[bel] = -0.5
        solDnh1[bel] = 0.5
    end
    for bel in extract_elements(select_boundaries(mesh, 1003))
        solDnh2[bel] = 0.5
    end
    for bel in extract_elements(select_boundaries(mesh, 1004))
        solDnf[bel] = 0.5
        solDnh1[bel] = -0.5
    end
    for bel in extract_elements(select_boundaries(mesh, 1005))
        solDnh3[bel] = -0.5
    end
    for bel in extract_elements(select_boundaries(mesh, 1006))
        solDnh3[bel] = 0.5
    end

    Dnf = compute_normalderivative(mesh, select_boundaries(mesh), 0.5 .* fh)
    any(abs.(Dnf[1] .- solDnf) .> 1e-15) && return false

    Dnh = compute_normalderivative(mesh, select_boundaries(mesh), hh, qdim=3)
    any(abs.(Dnh[1] .- solDnh1) .> 1e-15) && return false
    any(abs.(Dnh[2] .- solDnh2) .> 1e-15) && return false
    any(abs.(Dnh[3] .- solDnh3) .> 1e-15) && return false

    return true
end

@test test_utility()