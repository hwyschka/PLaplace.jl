"""
$(TYPEDSIGNATURES)
    
Dummy function for do nothing.
"""
function emptyfunction(args...) end

"""
    assemble_derivativetensor(
        mesh::Mesh;
        qdim::Int64 = 1
    ) -> Dict{Tuple{Int64, Int64}, SparseArrays.SparseMatrixCSC{Float64, Int64}}
    
Returns the discrete derivative tensor for all elements of mesh
and the number of components of the function qdim.
Each key tuple (j,r) represents the derivative matrix for the component r
in direction x_j.
"""
function assemble_derivativetensor(
    mesh::Mesh;
    qdim::Int64 = 1
)
    D = Dict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}()
    
    for j = 1:mesh.d
        for r = 1:qdim
            AA = zeros(Float64, mesh.nelems * (mesh.d+1))
            II = zeros(Int64, length(AA))
            JJ = zeros(Int64, length(AA))
            n = 0

            for el = 1:mesh.nelems
                nodes = mesh.Elements[el]
                (_, J) = jacobian(mesh, nodes)

                for k = 1:(mesh.d+1)
                    n = n + 1
                    II[n] = el
                    JJ[n] = qdim * (nodes[k]-1) + r
                    AA[n] = (J*grad_phi(mesh.d,k))[j]
                end
            end
            D[(j,r)] = sparse(II[1:n], JJ[1:n], AA[1:n], mesh.nelems, mesh.nnodes*qdim)
            dropzeros!(D[(j,r)])
        end
    end
    return D
end

"""
    assemble_derivativetensor_boundary(
        mesh::Mesh,
        BoundaryElements::Set{Int64}; 
        qdim::Int64 = 1
    ) -> Dict{Tuple{Int64, Int64}, SparseArrays.SparseMatrixCSC{Float64, Int64}}
    
Returns the discrete derivative tensor for all or specified boundary elements of mesh
and number of components qdim. Each key tuple (j,r) represents the derivative matrix for the
component r in direction x_j. Workaround based on the derivate tensor for the
corresponding full element and the fact that the gradient is constant on the element. 
"""
function assemble_derivativetensor_boundary(
    mesh::Mesh,
    BoundaryElements::Set{Int64}; 
    qdim::Int64 = 1
)
    D = Dict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}()
    
    for j = 1:mesh.d
        for r = 1:qdim
            AA = zeros(Float64, qdim * mesh.nboundelems * (mesh.d+1))
            II = zeros(Int64, length(AA))
            JJ = zeros(Int64, length(AA))
            n = 0

            for bel in BoundaryElements
                nodes = mesh.Elements[mesh.ParentElements[bel]]
                (_, J) = jacobian(mesh, nodes)

                for k = 1:(mesh.d+1)
                    n = n + 1
                    II[n] = bel
                    JJ[n] = qdim * (nodes[k]-1) + r
                    AA[n] = (J*grad_phi(mesh.d,k))[j]
                end
            end
            D[(j,r)] = sparse(II[1:n], JJ[1:n], AA[1:n], mesh.nboundelems, mesh.nnodes*qdim)
            dropzeros!(D[(j,r)])
        end
    end
    return D
end

"""
    assemble_derviativetensor_modified(
        D::AbstractDict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}},
        nodesToDrop::Set{Int64};
        qdim::Int64 = 1
    ) -> Dict{Tuple{Int64, Int64}, SparseArrays.SparseMatrixCSC{Float64, Int64}}
    
Returns a new derivative tensor reduced by nodesToDrop,
e.g used to drop homogeneous Dirichtlet boundary points of the mesh.
"""
function assemble_derviativetensor_modified(
    D::AbstractDict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}},
    nodesToDrop::Set{Int64};
    qdim::Int64 = 1
)
    DI = Dict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}()
    columnsToDrop = Set{Int64}()

    if qdim == 1
        columnsToDrop = nodesToDrop
    else
        for node in nodesToDrop
            union!(columnsToDrop, qdim*(node-1)+1:qdim*node)
        end
    end
    
    for (key, val) in D
        DI[key] = val[1:end, 1:end .∉ [columnsToDrop]]
    end
    return DI
end

"""
    compute_derivative(
        D::AbstractDict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}},
        v::AbstractVector{Float64}
    ) -> Dict{Tuple{Int64, Int64}, AbstractVector{Float64}}

Returs coefficient vectors of all the first partial derivatives ``D^{(j,r)}v`` of a function
given as coefficient vector v using the discrete derivative tensor D.
Each key tuple (j,r) represents the derivative of the component r in direction x_j.
"""
function compute_derivative(
    D::AbstractDict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}},
    v::AbstractVector{Float64}
)
    dv = Dict{Tuple{Int64,Int64},AbstractVector{Float64}}()
    for (key, val) in D
        dv[key] = val * v
    end
    return dv
end

"""
    compute_normalderivative(
        mesh::Mesh,
        boundary::Set{Boundary},
        v::AbstractVector{Float64};
        qdim::Int64 = 1
    ) -> Dict{Tuple{Int64, Int64}, AbstractVector{Float64}}

Returns normal derivative of a function on all or specified boundary elements of the mesh.
Each key r represents the derivative of the component r in outer normal direction.
"""
function compute_normalderivative(
    mesh::Mesh,
    BoundaryElements::Set{Int64},
    v::AbstractVector{Float64},
    n::AbstractVector{Float64};
    qdim::Int64 = 1
)
    Dnv = Dict{Int64,Vector{Float64}}()
    
    D = assemble_derivativetensor_boundary(mesh, BoundaryElements, qdim=qdim)
    dv = compute_derivative(D,v)

    nres = reshape(n, mesh.d, mesh.nboundelems)

    for r = 1:qdim
        Dnv[r] = zeros(Float64, mesh.nboundelems)
        for j = 1:mesh.d
            Dnv[r] += nres[j,:] .* dv[(j,r)]
        end
    end

    return Dnv
end

"""
    compute_normalderivative(
        mesh::Mesh,
        boundary::Set{Boundary},
        v::AbstractVector{Float64};
        qdim::Int64 = 1
    ) -> Dict{Tuple{Int64, Int64}, AbstractVector{Float64}}

Same as previous `$(FUNCTIONNAME)(...)`. However takes `Set{Boundary}` of named
boundary objects as argument `boundary` and translates it to a `Set{Int64}` containing
the indices of the nodes contained in the given boundaries.
"""
function compute_normalderivative(
    mesh::Mesh,
    boundary::Set{Boundary},
    v::AbstractVector{Float64};
    qdim::Int64 = 1
)
    boundary = extract_elements(boundary)
    return compute_normalderivative(
        mesh,
        boundary,
        v,
        outernormalvector(mesh, boundaryElements=boundary),
        qdim=qdim
    )
end


"""
$(TYPEDSIGNATURES)

Returns ``\\Vert f \\Vert^p_{X_p(\\Omega)} =
\\int_{\\Omega} \\Vert ∇ f(x) \\Vert_2^p \\;\\mathrm{d}x``
in FEM representaion, i.e.

```math
\\Vert f \\Vert^p_{X_p(T_{h_\\Omega})} =
    \\sum\\limits_{i=1}^m \\omega_i 
    \\left(\\sum\\limits_{j=1}^d \\sum\\limits_{r=1}^{d^\\prime}
    [D^{(j,r)} v]_i^2 \\right)^{\\frac{p}{2}}.
```

Parameters can either be given as a function and a mesh
or with the function as FEM coefficient vector and a derivative tensor.

# Mandatory Arguments
- `p::Float64`: Parameter for the norm.
- `dv::AbstractDict{Tuple{Int64,Int64},AbstractVector{Float64}}`: Discrete derivative of
    function v given at the quadrature nodes of the elements in ``T_{h_\\Omega}``.
- `w::AbstractVector{Float64}`: Vector of weights corresponding to the quadrature nodes.
"""
function xpnorm(
    p::Float64,
    dv::AbstractDict{Tuple{Int64,Int64},AbstractVector{Float64}},
    w::AbstractVector{Float64}
) :: Float64
    t = zeros(Float64, length(w))

    for (key, val) in dv
        t += val.^2
    end

    if isinf(p)
        return maximum(t.^(1/2))
    else
        return sum(w .* t.^(p/2))^(1/p)
    end
end

"""
$(TYPEDSIGNATURES)

Does the same as previous `$(FUNCTIONNAME)(...)`,
but takes coefficient vector and derivate tensor as arguments
instead of derivative vector.

Requires new mandatory arguments
- `v::AbstractVector{Float64}`: Function discretely evaluated at the nodes of a mesh.
- `D::AbstractDict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}`: Discrete derivative
    tensor corresponding to the same nodes as v.
which replace
- `dv::AbstractDict{Tuple{Int64,Int64},AbstractVector{Float64}}`
"""
function xpnorm(
    p::Float64,
    v::AbstractVector{Float64},
    D::AbstractDict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}},
    w::AbstractVector{Float64}
)
    dv = compute_derivative(D, v)
    return xpnorm(p, dv, w)
end

"""
    xpnorm(
        p::Float64,
        f::Function,
        mesh::Mesh;
        qdim::Int64 = 1
    ) -> Float64

Does the same as previous `$(FUNCTIONNAME)(...)`,
but takes a continous function and a mesh as arguments
instead of a discrete derivative (tensor) and coefficient vector.

Requires new mandatory arguments
- `f::Function`: Analytic description of the function to be evaluated.
- `mesh::Mesh`: FEM mesh corresponding to ``T_{h_\\Omega}``.
which replace
- `dv::AbstractDict{Tuple{Int64,Int64},AbstractVector{Float64}}`
"""
function xpnorm(
    p::Float64,
    f::Function,
    mesh::Mesh;
    qdim::Int64 = 1
)
    v = evaluate_mesh_function(mesh, f)
    D = assemble_derivativetensor(mesh, qdim=qdim)
    dv = compute_derivative(D, v)
    w = assemble_weightmultivector(mesh, qdim=qdim)
    return xpnorm(p, dv, w)
end

"""
    compute_lipschitzconstant_boundary(
        gh::Array{Float64,1},
        coords::Array{Array{Float64,1},1},
        boundaryNodes::Set{Int64};
        qdim::Int64 = 1,
        pnorm::Real = 2
    )

Returns Lipschitz constant of function g given as FEM coefficient vector
on a boundary of the mesh.

# Mandatory Arguments
- `gh::Array{Float64,1}`: Function discretely evaluated on a mesh.
- `coords::Array{Array{Float64,1},1}`: Array of coordinates of the nodes gh was evaluated
    on. In particular contains the coordinates of the boundaryNodes.
- `boundaryNodes::Set{Int64}`: Set of nodes, i.e. entries in coords that belong to the
    boundary that shall be evaluated. Therefore also selects entries of gh. Potentially
    shifted by qdim.

# Keyword Arguments
- `qdim::Int64 = 1`: Number of components of the function gh.
- `pnorm::Real = 2`: Parameter of p-norm used to measure internal distances.
"""
function compute_lipschitzconstant_boundary(
    gh::Array{Float64,1},
    coords::Array{Array{Float64,1},1},
    boundaryNodes::Set{Int64};
    qdim::Int64 = 1,
    pnorm::Real = 2
)
    sup::Float64 = 0
    for x in boundaryNodes
        for y in boundaryNodes
            var = norm(gh[(qdim*(y-1)+1):qdim*y] - gh[(qdim*(x-1)+1):qdim*x], pnorm) / 
                norm(coords[y] .- coords[x], pnorm)
            if var > sup
                sup = var
            end
        end
    end

    return sup
end
