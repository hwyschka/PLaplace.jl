"""
$(TYPEDEF)

Structure to store statistics to the algorithm during the runtime.
The information is supposed to be exported to [PLaplaceData](@ref) for external usage.

# Fields
$(TYPEDFIELDS)
"""
mutable struct AlgorithmData   
    "Obtained accuracy."
    eps::Float64

    "Required iterations for the auxilliary path-following."
    Naux::Union{Int64,Missing}

    "Required iterations for the main path-following."
    Nmain::Union{Int64,Missing}

    "Required time for the setup."
    tsetup::Union{Float64,Missing}

    "Required time for the auxilliary path-following."
    taux::Union{Float64,Missing}

    "Required time for the main path-following."
    tmain::Union{Float64,Missing}

    "Solution vector if iteration converged."
    solution::Union{Vector{Float64},Missing}

    "Notifications from the iteration."
    msg::String
end


"""
$(TYPEDSIGNATURES)

Constructor for [AlgorithmData](@ref) with empty values.
"""
function AlgorithmData()
    return AlgorithmData(Inf, missing, missing, missing, missing, missing, missing, "-")
end

"""
$(TYPEDSIGNATURES)

Handling the termination of the algorithm because the maximum number of iterations during
a phase was reached.
"""
function handle_maxiterations!(data::AlgorithmData, phase::String, iteration::Int64)
    data.msg *= " Exceeded iterations in $phase$iteration -"
end

"""
$(TYPEDSIGNATURES)

Handling the termination of the algorithm because the stepsize update parameter κ in a 
path-following with adaptive stepping got numerically too small.
Should theoretically not occur, so this usually indicates an infeasible problem.
"""
function handle_kappavanish!(data::AlgorithmData, phase::String, iteration::Int64)
    data.msg *= " ⲕ too small in $phase$iteration -"
end

"""
$(TYPEDSIGNATURES)

Handling the termination of the algorithm because the final update step in the auxilliary
path-following failed because of a singular system matrix. 
"""
function handle_auxfail!(data::AlgorithmData)
    data.msg *= " Auxiliary criteria not fulfilled  -"
end

"""
$(TYPEDSIGNATURES)

Handling the termination of a main path-following because of a singular system matrix. 
In particular stores the resulting (reversly computed) obtained accuracy.
"""
function handle_accuracy!(
    data::AlgorithmData,
    eps::Float64
)
    data.eps = eps
end

"""
$(TYPEDSIGNATURES)

Handling the assembly of barrier terms. 
In particular stores message if solver or preconditioner got changed.
"""
function handle_assembly!(
    data::AlgorithmData,
    log::LogData,
    tracker::AssemblyTracker,
    phase::String,
    iteration::Int64
)
    if tracker.factorization
        data.msg *= " Changed to LU fact. in $phase$iteration -"
        log_change_factorization(log)
    end

    if tracker.singularity
        data.msg *= " Hessian singular in $phase$iteration -"
    end

    if tracker.preconditioner
        data.msg *= " Changed to LU prec. in $phase$iteration -"
        log_change_preconditioner(log)
    end
end

"""
$(TYPEDSIGNATURES)
    
Wrapper for passing statistics to the log handling.
"""
function log_statistics(A::AlgorithmData, L::LogData)
    return log_statistics(
        L,
        A.tsetup,
        A.taux,
        A.tmain,
        A.Naux,
        A.Nmain,
        A.eps,
        A.msg
    )
end
