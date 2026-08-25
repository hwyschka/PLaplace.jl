"""
$(TYPEDEF)

Object used to track the computation of the dual norm described
by the Hessian of a self-concordant barrier.

# Fields
$(TYPEDFIELDS)
"""
mutable struct StarnormTracker  
    "Flag when Hessian was indefinite and result of the scalar product negative"
    failed::Bool
end


"""
$(TYPEDSIGNATURES)

Constructor for creating a default [StarnormTracker](@ref) object.
"""
function StarnormTracker() 
    return StarnormTracker(false)
end

"""
$(TYPEDSIGNATURES)

Sets iteration count, final scaling value and failed flag to a [StarnormTracker](@ref)
when no admissible descent direction could be obtained in the given iteration limit.
"""
function fail!(tracker::StarnormTracker)
    tracker.failed = true
end

"""
$(TYPEDSIGNATURES)

Resets a [StarnormTracker](@ref) to its default constructor values.
"""
function reset!(tracker::StarnormTracker)
    tracker.failed = false
end

"""
$(TYPEDEF)

Object used to track the result of the application of a descent step.
Usually done by a backtracking line search, it contains the required
iterations, i.e. updates of the scaling parameter, the final scaling
parameter and the required time.

# Fields
$(TYPEDFIELDS)
"""
mutable struct DescentTracker
    "Number of iterations."
    i::Int64
    
    "Final scaling value"
    val::Float64

    "Flag descent failed for one of the tracked reasons"
    failed::Bool
    
    "Flag when starnorm for damping could not be obtained"
    nodamping::Bool

    "Flag when no admissible direction could be obtained"
    nodescent::Bool

    "Required time"
    time::Union{Float64,Missing}
end


"""
$(TYPEDSIGNATURES)

Constructor for creating a default [DescentTracker](@ref) object.
"""
function DescentTracker() 
    return DescentTracker(0, 1.0, false, false, false, missing)
end

"""
$(TYPEDSIGNATURES)

Constructor for creating a [DescentTracker](@ref) object
with given iteration count and scaling parameter.
"""
function DescentTracker(i::Int64, val::Float64)
    return DescentTracker(i, val, false, false, false, missing)
end

"""
$(TYPEDSIGNATURES)

Sets iteration count and final scaling value to a [DescentTracker](@ref).
"""
function set!(tracker::DescentTracker, i::Int64, val::Float64)
    tracker.i = i
    tracker.val = val
end

"""
$(TYPEDSIGNATURES)

Sets iteration count, final scaling value and failed flag to a [DescentTracker](@ref)
when no admissible descent direction could be obtained in the given iteration limit.
"""
function fail_direction!(tracker::DescentTracker, i::Int64, val::Float64)
    tracker.i = i
    tracker.val = val
    tracker.nodescent = true
    tracker.failed = true
end

"""
$(TYPEDSIGNATURES)

Sets failed flag to a [DescentTracker](@ref)
when the dual norm for the damping parameter could not be computed.
"""
function fail_damping!(tracker::DescentTracker)
    tracker.i = 0
    tracker.val = 1.0
    tracker.nodamping = true
    tracker.failed = true
end

"""
$(TYPEDSIGNATURES)

Resets a [DescentTracker](@ref) to its default constructor values.
"""
function reset!(tracker::DescentTracker)
    tracker.i = 0
    tracker.val = 1.0
    tracker.nodamping = false
    tracker.nodescent = false
    tracker.failed = false
    tracker.time = missing
end

"""
$(TYPEDEF)
    
Tracks if factorization or preconditioner got changed
during an assembly of the IterationData.

Also tracks condition of system of system matrix, i.e. the barrier Hessian.
Has to be tracked via this tool because after the assembly
the matrix is usually only available as factorization.

For debugging purposes also the explicit hessian can be exported,
but this becomes slow for larger problems. 

# Fields
$(TYPEDFIELDS)
"""
mutable struct AssemblyTracker
    "Flag if conditino is tracked."
    trackcondition::Bool

    "Condition number of system matrix."
    conditionnumber::Union{Float64,Missing}

    "Flag if full Hessian is tracked."
    trackhessian::Bool

    "Full system matrix"
    hessian::Union{SparseMatrixCSC{Float64, Int64},Missing}

    "Flag if hessian is singular."
    singularity::Bool
    
    "Flag if factorization got changed."
    factorization::Bool
    
    "Flag if precondtioner got changed."
    preconditioner::Bool
end

"""
    AssemblyTracker(;trackcondition::Bool = false)

Default constructor for [AssemblyTracker](@ref).
"""
function AssemblyTracker(;trackcondition::Bool = false, trackhessian::Bool = false)
    return AssemblyTracker(
        trackcondition,
        missing,
        trackhessian,
        missing,
        false,
        false,
        false
    )
end

"""
$(TYPEDSIGNATURES)

Resets an [AssemblyTracker](@ref) to its default constructor values.
Does not change outside flag if condition number is tracked.
"""
function reset!(tracker::AssemblyTracker)
    tracker.conditionnumber = missing
    tracker.hessian = missing
    tracker.singularity = false
    tracker.factorization = false
    tracker.preconditioner = false
end

"""
$(TYPEDSIGNATURES)

Sets flag of an [AssemblyTracker](@ref) that the hessian is singular
and no factorization was computed.
"""
function hessian_singular!(tracker::AssemblyTracker)
    tracker.singularity = true
end

"""
$(TYPEDSIGNATURES)

Sets flag of an [AssemblyTracker](@ref) that the factorization got changed.
"""
function changed_factorization!(tracker::AssemblyTracker)
    tracker.factorization = true
end

"""
$(TYPEDSIGNATURES)

Sets flag of an [AssemblyTracker](@ref) that the preconditioner got changed.
"""
function changed_preconditioner!(tracker::AssemblyTracker)
    tracker.preconditioner = true
end
