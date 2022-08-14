"""
    emptyfunction(args...)
    
Dummy function for do nothing.
"""
function emptyfunction(args...) end

"""
    assemble_derivativetensor(mesh::Mesh) -> Dict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}
    
Returns the discrete derivative tensor for all elements of mesh and image dimension qdim.
Each key tuple (j, r) represents the derivative matrix for the image dimension r in direction x_j.
"""
function assemble_derivativetensor(mesh::Mesh; qdim=1)
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
    assemble_derivativetensor_boundary(mesh::Mesh, BoundaryElements::Set{Int64}; qdim=1) -> Dict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}
    
Returns the discrete derivative tensor for all or specied boundary elements of mesh and image dimension qdim.
Each key tuple (j, r) represents the derivative matrix for the image dimension r in direction x_j.
Workaround based on the derivate tensor for the corresponding full element
and the fact that the gradient ist constant on the element. 
"""
function assemble_derivativetensor_boundary(mesh::Mesh, BoundaryElements::Set{Int64}; qdim=1)
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
    assemble_derviativetensor_modified(D::AbstractDict{Int64,SparseMatrixCSC{Float64,Int64}}, nodesToDrop::Set{Int64}; qdim=1) -> Dict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}
    
Returns a new derivative tensor reduced by nodesToDrop, e.g used to drop Dirichtlet 0 boundary points of the mesh.
"""
function assemble_derviativetensor_modified(D::AbstractDict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}, nodesToDrop::Set{Int64}; qdim=1)
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
    compute_derivative(D::AbstractDict{Int64,SparseMatrixCSC{Float64,Int64}}, v::AbstractVector{Float64}) -> Dict{Tuple{Int64,Int64},AbstractVector{Float64}}()

Returs coefficient vectors of the first partial derivatives Dv[j] of a function given 
as coefficient vector v using the discrete derivative tensor D.
Each key tuple (j, r) represents the derivative of the image dimension r in direction x_j.
"""
function compute_derivative(D::AbstractDict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}, v::AbstractVector{Float64})
    dv = Dict{Tuple{Int64,Int64},AbstractVector{Float64}}()
    for (key, val) in D
        dv[key] = val * v
    end
    return dv
end

"""
    compute_normalderivative(mesh::Mesh, boundary::Set{Boundary}, v::AbstractVector{Float64}; qdim::Int64=1) -> Dict{Int64,Vector{Float64}}()
    compute_normalderivative(mesh::Mesh, BoundaryElements::Set{Int64}, v::AbstractVector{Float64}, n::AbstractVector{Float64}; qdim::Int64=1) -> Dict{Int64,Vector{Float64}}()

Returns normal derivative of a function on all or specified boundary elements of the mesh.
Each key r represents the derivative of the image dimension r in outer normal direction.
"""
function compute_normalderivative(mesh::Mesh, boundary::Set{Boundary}, v::AbstractVector{Float64}; qdim::Int64=1)
    boundary = extract_elements(boundary)
    return compute_normalderivative(mesh, boundary, v, outernormalvector(mesh, boundaryElements=boundary), qdim=qdim)
end

function compute_normalderivative(mesh::Mesh, BoundaryElements::Set{Int64}, v::AbstractVector{Float64}, n::AbstractVector{Float64}; qdim::Int64=1)
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
    xpnorm(p::Float64, dv::AbstractDict{Int64,AbstractVector{Float64}}, w::AbstractVector{Float64}, m::Int64, dim::Int8, qdim::Int64) -> Vector{Float64}
    xpnorm(p::Float64, v::AbstractVector{Float64}, D::AbstractDict{Int64,SparseMatrixCSC{Float64,Int64}}, w::AbstractVector{Float64}, m::Int64, qdim::Int64) -> Vector{Float64}
    xpnorm(p::Float64,f::Function, mesh::Mesh; qdim=1) -> Vector{Float64}

Returns ``||f||^p_{X_p} = int_{Omega} norm{∇ f(x)}_2^p``
in FEM representaion, i.e. 
``sum_{i=1}^m omega_i (sum_{j=1}^d ((D^{(j)} v)_i)^2)^{p/2}``.

Parameters can either be given as a function and a mesh
or with the function as FEM coefficient vector and a derivative tensor.

"""
function xpnorm(p::Float64, dv::AbstractDict{Tuple{Int64,Int64},AbstractVector{Float64}}, w::AbstractVector{Float64})
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

function xpnorm(p::Float64, v::AbstractVector{Float64}, D::AbstractDict{Tuple{Int64,Int64},SparseMatrixCSC{Float64,Int64}}, w::AbstractVector{Float64})
    dv = compute_derivative(D, v)
    return xpnorm(p, dv, w)
end

function xpnorm(p::Float64,f::Function, mesh::Mesh; qdim=1)
    v = evaluate_mesh_function(mesh, f)
    D = assemble_derivativetensor(mesh, qdim=qdim)
    dv = compute_derivative(D, v)
    w = assemble_weightmultivector(mesh, qdim=qdim)
    return xpnorm(p, dv, w)
end
