"""
$(TYPEDSIGNATURES)

Returns discrete source term based on the volume source `f`.
`f` has `qdim` components and can be given as analytical function or
discrete coefficient vector on the nodes or quadrature points of `mesh`.
"""
function volume_source(
    mesh::Mesh,
    f::Union{AbstractVector{Float64}, Function, Missing},
    qdim::Int64
) :: Vector{Float64}
    src = zeros(Float64, qdim*mesh.nnodes)

    if !ismissing(f) && !iszero(f)
        if f isa Function
            fd =  evaluate_mesh_function(mesh, f, qdim=qdim)
        else
            fd = f            
        end

        if length(fd) == mesh.nnodes*qdim
            M = assemble_massmatrix(
                mesh,
                qdim = qdim,
                order = 3
            )
            src += M * fd
        elseif mod(length(fd), mesh.nelems*qdim) == 0
            nPoints = div(length(fd), mesh.nelems * qdim)
            quadOrder = quadrature_order(mesh.d, nPoints)
            
            E = assemble_basismatrix(
                mesh,
                qdim = qdim,
                order = quadOrder
            )
            W = Diagonal(
                assemble_weightmultivector(
                    mesh,
                    qdim = qdim,
                    order = quadOrder
                )
            )
            src += E' * W * fd
        else
            throw(DomainError(fd,"Dimension Missmatch"))
        end
    end

    return src
end

"""
$(TYPEDSIGNATURES)

Returns discrete source term based on the boundary source `h`.
`h` has `qdim` components and can be given as analytical function or
discrete coefficient vector on the nodes or quadrature points of `mesh`.
Note that if `h` is given as an analytical function, it will be evaluated at the mesh nodes
where properties like the normal vector might not be defined.
If you want to account for such potential issues,
you have to manually evaluate at quadrature points before.
"""
function boundary_source(
    mesh::Mesh,
    neumann_boundary::Union{Set{Boundary}, Set{Int64}, Missing},
    h::Union{AbstractVector{Float64}, Function, Missing},
    qdim::Int64
) :: Vector{Float64}
    src = zeros(Float64, qdim*mesh.nnodes)

    if !ismissing(neumann_boundary) && !ismissing(h) && !iszero(h)
        if neumann_boundary isa Set{Boundary}
            belems = extract_elements(neumann_boundary)
            bnodes = extract_nodes(neumann_boundary)
        else
            belems = neumann_boundary
            bnodes = Set{Int64}()
            if !ismissing(belems)
                for el in belems
                    for node in mesh.BoundaryElements[el]
                        push!(bnodes, node)
                    end
                end
            end
        end

        if h isa Function
            hd = evaluate_mesh_function(mesh, h, region=_neumann_nodes, qdim=qdim)
        else
            hd = h            
        end
        
        if length(hd) == mesh.nnodes*qdim
            N = assemble_massmatrix_boundary(
                mesh,
                boundaryElements = belems,
                qdim = qdim,
                order = 3
            )
            src += N * hd
        elseif mod(length(hd), mesh.nboundelems * qdim) == 0
            nPoints = div(length(hd), mesh.nboundelems * qdim)
            quadOrder = quadrature_order(mesh.d - 1, nPoints)

            E = assemble_basismatrix_boundary(
                mesh,
                boundaryElements = belems,
                qdim = qdim,
                order = quadOrder
            )
            W = Diagonal(
                assemble_weightmultivector_boundary(
                    mesh,
                    qdim = qdim,
                    order = quadOrder
                )
            )
            src += E' * W * hd
        else
            throw(DomainError(hd,"Dimension Missmatch"))
        end
    end

    return src
end

"""
$(TYPEDSIGNATURES)

Returns value of source term based on the volume source `f` and the boundary source `h`.
For more details see [volume_source](@ref) and [boundary_source](@ref).
"""
function sources_terms(
    mesh::Mesh,
    neumann_boundary::Union{Set{Boundary}, Set{Int64}, Missing},
    h::Union{AbstractVector{Float64}, Function, Missing},
    f::Union{AbstractVector{Float64}, Function, Missing},
    qdim::Int64
) :: Vector{Float64}
    src = zeros(Float64, qdim*mesh.nnodes)
    src += volume_source(mesh, f, qdim)
    src += boundary_source(mesh, neumann_boundary, h, qdim)

    return src
end

"""
$(TYPEDSIGNATURES)

Returns value of source term based on the volume source `f` and the boundary source `h`
evaluated at `u`.
For more details see [sources_terms](@ref).
"""
function compute_sources(
    u::AbstractVector{Float64},
    mesh::Mesh,
    neumann_boundary::Union{Set{Boundary}, Set{Int64}, Missing},
    h::Union{AbstractVector{Float64}, Function, Missing},
    f::Union{AbstractVector{Float64}, Function, Missing},
    qdim::Int64
) :: Float64
    s = sources_terms(mesh, neumann_boundary, h, f, qdim)
    return dot(s,u)
end

"""
$(TYPEDSIGNATURES)

Computes value of characteristic derivative term in p-Laplace functional evaluated at `u`.
"""
function compute_plaplace_term(
    u::AbstractVector{Float64},
    p::Float64,
    mesh::Mesh,
    qdim::Int64
) :: Float64
    D = assemble_derivativetensor(mesh, qdim=qdim)

    Du = zeros(Float64, mesh.nelems)
    for (key,val) in D
        Du += (val*u).^2
    end

    w = assemble_weightmultivector(mesh, qdim=1, order=1)
    return (1.0 / p) * dot(w, Du.^(p/2)) 
end

"""
    objective_functional(
        u::AbstractVector{Float64},
        p::Float64,
        mesh::Mesh; 
        f::Union{AbstractVector{Float64}, Function, Missing} = missing,
        neumann_boundary::Union{Set{Boundary}, Set{Int64}, Missing} = missing,
        h::Union{AbstractVector{Float64}, Function, Missing} = missing,
        qdim::Int64 = 1
    ) -> Float64

Returns value of variational formulation functional for the p-Laplace problem
evaluated at `u`. 
"""
function objective_functional(
    u::AbstractVector{Float64},
    p::Float64,
    mesh::Mesh; 
    f::Union{AbstractVector{Float64}, Function, Missing} = missing,
    neumann_boundary::Union{Set{Boundary}, Set{Int64}, Missing} = missing,
    h::Union{AbstractVector{Float64}, Function, Missing} = missing,
    qdim::Int64 = 1
) :: Float64
    if p < 1
        throw(DomainError(p, "This package only supports 1 ≤ p ≤ ∞."))
    end

    s = compute_sources(
        u,
        mesh,
        neumann_boundary,
        h,
        f,
        qdim
    )

    t = compute_plaplace_term(u, p, mesh, qdim)

    return t - s
end

"""
$(TYPEDSIGNATURES)

Returns various errors of numerical solution compared to given 
discretized analytical solution. 
First value is error in difference in objective value,
then pointwise error L1, L2 and LInf norm.

Intended to be called via a wrapper to ensure the discretized analytical solution
and the numerical solution match the mesh.
"""
function compute_errors(
    anasol::Union{AbstractVector{Float64}, Function},
    numsol::AbstractVector{Float64},
    mesh::Mesh,
    neumann_boundary::Union{Set{Boundary}, Set{Int64}, Missing},
    h::Union{AbstractVector{Float64}, Missing},
    f::Union{AbstractVector{Float64}, Missing},
    p::Float64,
    qdim::Int64
) :: Array{Float64,1}
    error = anasol .- numsol 

    o1 = objective_functional(
        anasol,
        p,
        mesh,
        f = f,
        neumann_boundary = neumann_boundary,
        h = h,
        qdim = qdim
    )
    o2 = objective_functional(
        numsol,
        p,
        mesh,
        f = f,
        neumann_boundary = neumann_boundary,
        h = h,
        qdim = qdim
    )

    errors = [
        o2-o1,
        pnorm(1.0, error, mesh, qdim=qdim, order=3),
        pnorm(2.0, error, mesh, qdim=qdim, order=3),
        pnorm(Inf, error, mesh, qdim=qdim, order=3)
    ]

    return errors
end
