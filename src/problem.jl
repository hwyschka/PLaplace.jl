"""
    objective_functional(
        u::AbstractVector{Float64},
        p::Float64,
        mesh::Mesh; 
        f::AbstractVector{Float64}=Vector{Float64}(),
        neumannBoundary::Set{Boundary}=Set{Boundary}(),
        h::AbstractVector{Float64}=Vector{Float64}(),
        qdim::Int64=1
    )
    objective_functional(
        u::Function,
        p::Float64,
        mesh::Mesh; 
        f::Function=emptyfunction,
        neumannBoundary::Set{Boundary}=Set{Boundary}(),
        h::Function=emptyfunction,
        qdim::Int64=1
    )

Returns value of variational formulation functional for the p-Laplace problem
evaluated at u. 
"""
function objective_functional(
    u::AbstractVector{Float64},
    p::Float64,
    mesh::Mesh; 
    f::AbstractVector{Float64}=Vector{Float64}(),
    neumannBoundary::Set{Boundary}=Set{Boundary}(),
    h::AbstractVector{Float64}=Vector{Float64}(),
    qdim::Int64=1
)
    if p < 1
        throw(DomainError(p, "This package only supports 1 ≤ p ≤ ∞."))
    end

    s = compute_sources(
        u,
        mesh,
        neumannBoundary,
        h,
        f,
        qdim
    )

    t = compute_plaplace_term(u, p, mesh, qdim)

    return t + s
end

function objective_functional(
    u::Function,
    p::Float64,
    mesh::Mesh; 
    f::Function=emptyfunction,
    neumannBoundary::Set{Boundary}=Set{Boundary}(),
    h::Function=emptyfunction,
    qdim::Int64=1
)
    _u = evaluate_mesh_function(mesh, u, qdim=qdim)
    _f = f == emptyfunction ? 
        zeros(Float64, mesh.nnodes * qdim) :
        evaluate_mesh_function(mesh, f, qdim=qdim) 
    _h = (isempty(neumannBoundary) || h == emptyfunction) ?
        zeros(Float64, mesh.nnodes * qdim) :
        evaluate_mesh_function(mesh, h, neumannBoundary, qdim=qdim)

    return objective_functional(
        _u,
        p,
        mesh,
        f=_f,
        neumannBoundary=neumannBoundary,
        h=_h, 
        qdim=qdim
    )
end

"""
    compute_plaplace_term(
        u::AbstractVector{Float64},
        p::Float64,
        mesh::Mesh,
        qdim::Int64
    )

Computes value of characteristic derivative term in p-Laplace functional evaluated at u.
"""
function compute_plaplace_term(
    u::AbstractVector{Float64},
    p::Float64,
    mesh::Mesh,
    qdim::Int64
)
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
    compute_sources(
        u::AbstractVector{Float64},
        mesh::Mesh,
        neumannBoundary::Set{Boundary},
        h::AbstractVector{Float64},
        f::AbstractVector{Float64},
        qdim::Int64
    )    

Computes value of source terms in p-Laplace functional evaluated at u.
"""
function compute_sources(
    u::AbstractVector{Float64},
    mesh::Mesh,
    neumannBoundary::Set{Boundary},
    h::AbstractVector{Float64},
    f::AbstractVector{Float64},
    qdim::Int64
)    
    rhs = zeros(Float64, qdim*mesh.nnodes)
    
    if any(o -> o != 0, f)
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
    
    if !isempty(neumannBoundary)
        if length(h) == mesh.nnodes*qdim
            N = assemble_massmatrix_boundary(
                mesh,
                boundaryElements=extract_elements(neumannBoundary),
                qdim=qdim,
                order=3
            )
            rhs -= N * h
        elseif mod(length(h),mesh.nboundelems*qdim) == 0
            nPoints = div(length(h), mesh.nboundelems * qdim)
            quadOrder = quadrature_order(mesh.d-1, nPoints)

            E = assemble_basismatrix_boundary(
                mesh,
                boundaryElements=extract_elements(neumannBoundary),
                qdim=qdim,
                order=quadOrder
            )
            W = Diagonal(
                assemble_weightmultivector_boundary(
                    mesh,
                    qdim=qdim,
                    order=quadOrder
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
    compute_errors(
        anasol::Array{Float64,1},
        numsol::Array{Float64,1},
        mesh::Mesh,
        neumannBoundary::Set{Boundary},
        h::AbstractVector{Float64},
        f::AbstractVector{Float64},
        p::Float64,
        qdim::Int64
    )
    compute_errors(
        anasol::Array{Float64,1},
        outputData::PLaplaceData
    )

Returns various errors of numerical solution compared to given 
discretized analytical solution. 
First value is error in difference in objective value,
then pointwise error L1, L2 and LInf norm.
"""
function compute_errors(
    anasol::Array{Float64,1},
    numsol::Array{Float64,1},
    mesh::Mesh,
    neumannBoundary::Set{Boundary},
    h::AbstractVector{Float64},
    f::AbstractVector{Float64},
    p::Float64,
    qdim::Int64
)
    error = anasol .- numsol 

    o1 = objective_functional(
        anasol,
        p,
        mesh,
        f=f,
        neumannBoundary=neumannBoundary,
        h=h,
        qdim=qdim)
    o2 = objective_functional(
        numsol,
        p,
        mesh,
        f=f,
        neumannBoundary=neumannBoundary,
        h=h,
        qdim=qdim
    )

    errors = [
        o2-o1,
        pnorm(1.0, error, mesh, qdim=qdim, order=3),
        pnorm(2.0, error, mesh, qdim=qdim, order=3),
        pnorm(Inf, error, mesh, qdim=qdim, order=3)
    ]

    return errors
end

function compute_errors(
    anasol::Array{Float64,1},
    outputData::PLaplaceData
)
    return compute_errors(
        anasol,
        outputData.v,
        outputData.mesh,
        outputData.neumannBoundary,
        outputData.h,
        outputData.f,
        outputData.p,
        outputData.qdim
    )
end
