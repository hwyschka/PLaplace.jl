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

    "Required iterations for the auxiliary path-following."
    Naux::Union{Int64,Missing}

    "Required iterations for the main path-following."
    Nmain::Union{Int64,Missing}

    "Required time for the setup."
    tsetup::Union{Float64,Missing}

    "Required time for the auxiliary path-following."
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
function handle_maxiterations!(
    data::AlgorithmData,
    phase::String,
    iteration::Int64
)
    data.msg *= " Exceeded iterations in $phase$iteration -"
end

"""
$(TYPEDSIGNATURES)

Handling the termination of the algorithm because the stepsize update parameter κ in a 
path-following with adaptive stepping got numerically too small.
Should theoretically not occur, so this usually indicates an infeasible problem.
"""
function handle_kappavanish!(
    data::AlgorithmData,
    phase::String,
    iteration::Int64
)
    data.msg *= " ⲕ too small in $phase$iteration -"
end

"""
$(TYPEDSIGNATURES)

Handling the termination of the algorithm because the stepsize update parameter κ in a 
path-following with adaptive stepping got numerically too small.
Should theoretically not occur, so this usually indicates an infeasible problem.
"""
function handle_starnorm!(
    data::AlgorithmData,
    computation::String,
    phase::String,
    iteration::Int64
)
    data.msg *= " Norm indefinite for $computation in $phase$iteration -"
end

"""
$(TYPEDSIGNATURES)

Handling the termination of the algorithm because the final update step in the auxiliary
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
    if tracker.trackhessian
        log_debug_hessian(log, iteration, tracker.hessian)
    end
    
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

Handling the application of the descent direction.
In particular stores message if iterations in backtracking were exceeded.
Otherwise processes the assembly.
"""
function handle_descent!(
    data::AlgorithmData,
    log::LogData,
    tracker::Union{Missing, DescentTracker},
    assemblytracker::AssemblyTracker,
    phase::String,
    iteration::Int64
)
    if !ismissing(tracker) && tracker.nodamping
        handle_starnorm!(data, "damping", phase, iteration)
    elseif !ismissing(tracker) && tracker.nodescent
        data.msg *= " Backtracking failed in $phase$iteration -"
    else
        handle_assembly!(data, log, assemblytracker, phase, iteration)
    end
end
