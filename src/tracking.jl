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
    
    "Required time"
    time::Union{Float64,Missing}
end


"""
$(TYPEDSIGNATURES)

Constructor for creating a default [DescentTracker](@ref) object.
"""
function DescentTracker() 
    return DescentTracker(0, 0.0, missing)
end

"""
$(TYPEDSIGNATURES)

Constructor for creating a [DescentTracker](@ref) object
with given iteration count and scaling parameter.
"""
function DescentTracker(i::Int64, val::Float64)
    return DescentTracker(i, val, missing)
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

Resets a [DescentTracker](@ref) to its default constructor values.
"""
function reset!(tracker::DescentTracker)
    tracker.i = 0
    tracker.val = 0.0
    tracker.time = missing
end

"""
$(TYPEDEF)
    
Tracks if factorization or preconditioner got changed
during an assembly of the IterationData.

Also tracks condition of system of system matrix, i.e. the barrier Hessian.
Has to be tracked via this tool because after the assembly
the matrix is usually only available as factorization.

# Fields
$(TYPEDFIELDS)
"""
mutable struct AssemblyTracker
    "Flag if conditino is tracked."
    trackcondition::Bool

    "Condition number of system matrix."
    conditionnumber::Union{Float64,Missing}

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
function AssemblyTracker(;trackcondition::Bool = false)
    return AssemblyTracker(trackcondition, missing, false, false, false)
end

"""
$(TYPEDSIGNATURES)

Resets an [AssemblyTracker](@ref) to its default constructor values.
Does not change outside flag if condition number is tracked.
"""
function reset!(tracker::AssemblyTracker)
    tracker.conditionnumber = missing
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
