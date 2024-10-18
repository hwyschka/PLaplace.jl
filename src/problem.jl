"""
$(TYPEDSIGNATURES)

Computes value of characteristic derivative term in p-Laplace functional evaluated at u.
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

    if isinf(p)
        return maximum(Du.^(1/2))
    else
        w = assemble_weightmultivector(mesh, qdim=1, order=1)
        return (1.0 / p) * dot(w, Du.^(p/2)) 
    end
    
end

"""
$(TYPEDSIGNATURES)

Computes value of the required boundary and volume source terms for the p-Laplace
functional evaluated at u.
"""
function compute_sources(
    u::AbstractVector{Float64},
    mesh::Mesh,
    neumann_boundary::Union{Set{Boundary}, Set{Int64}, Missing},
    h::Union{AbstractVector{Float64}, Function, Missing},
    f::Union{AbstractVector{Float64}, Function, Missing},
    qdim::Int64
) :: Float64
    rhs = zeros(Float64, qdim*mesh.nnodes)
    
    if !ismissing(f) && any(o -> o != 0, f)
        if length(f) == mesh.nnodes*qdim
            M = assemble_massmatrix(
                mesh,
                qdim=qdim,
                order=3
            )
            rhs -= M * f
        elseif mod(length(f), mesh.nelems*qdim) == 0
            nPoints = length(f) / (mesh.nelems * qdim)
            quadOrder = quadrature_order(mesh.d, nPoints)
            
            E = assemble_basismatrix(
                mesh,
                qdim=qdim,
                order=quadOrder
            )
            W = Diagonal(
                assemble_weightmultivector(
                    mesh,
                    qdim=qdim,
                    order=quadOrder
                )
            )
            rhs -= E' * W * f
        else
            throw(DomainError(f,"Dimension Missmatch"))
        end
    end
    
    if !ismissing(neumann_boundary)
        if neumann_boundary isa Set{Boundary}
            belems = extract_elements(neumann_boundary)
        else
            belems = neumann_boundary
        end
        
        if length(h) == mesh.nnodes*qdim
            N = assemble_massmatrix_boundary(
                mesh,
                boundaryElements = belems,
                qdim = qdim,
                order = 3
            )
            rhs -= N * h
        elseif mod(length(h),mesh.nboundelems*qdim) == 0
            nPoints = div(length(h), mesh.nboundelems * qdim)
            quadOrder = quadrature_order(mesh.d-1, nPoints)

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
            rhs -= E' * W * h
        else
            throw(DomainError(h,"Dimension Missmatch"))
        end
    end

    return dot(u,rhs)
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
evaluated at u. 
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

    return t + s
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
    anasol::Array{Float64,1},
    numsol::Array{Float64,1},
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
