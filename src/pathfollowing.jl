"""
$(TYPEDEF)

Public type used to specify the stepping scheme of the path-following schemes.
For more information see the [interior-point section](@ref path-following-theory)
in the documentation.

# Available Options
- `SHORT`:
    Regular update of the paramter t.
- `LONG`:
    Larger update of the parameter t if iterate fulfill approximate centering conditon.
    Results in line-search for application of update on the iterate x.
- `ADAPTIVE`:
    Same as the long stepping, but the factor for the larger update is adaptively changed
    depending on how many steps it required to fulfill the approximate centering condition.
"""
@enum Stepsize begin
    SHORT
    LONG
    ADAPTIVE
end

"""
$(TYPEDSIGNATURES)

Executes auxilliary path-following with adaptive stepsize.

The iteration is performed on the [IterationData](@ref), which will later store
the final iterate as well as all the corresponding barrier terms.
The number of required iterations and potential messages will be stored in
[AlgorithmData](@ref).
If [LogData](@ref) is verbose, data per iteration will be written to the output stream.
Further, if a file is provided, the log will be also exported to that file.
"""
function pathfollowing_auxilliary_adaptive!(
    I::IterationData,
    A::AlgorithmData,
    S::StaticData,
    L::LogData
)
    searchtracker = DescentTracker(0,0.0)
    assemblytracker = AssemblyTracker(trackcondition=L.trackcondition)
    log_inital(L)

    assemble!(I, S, tracker=assemblytracker)
    handle_assembly!(A, L, assemblytracker,"A", 0)
    
    t::Float64 = 1
    kappa::Float64 = S.kappa
    maxIter::Int64 = S.maxIter
    
    iterationCount::Int64 = 0
    G::AbstractVector{Float64} = -I.gradientF
    bound::Float64 = sqrt(S.beta) / (1 + sqrt(S.beta))

    lastAccept::Int64 = 0
    lastAcceptedt::Float64 = t
    lastAcceptedx::AbstractVector{Float64} = I.x 

    if assemblytracker.singularity
        iterationCount = -1
        maxIter = 0
    else
        critnorm = starnorm(G, I, S.solveLS)

        log_iteration(
            L,
            0, 0 , "A",
            missing, missing, missing, missing,
            missing, missing,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if critnorm <= bound
            iterationCount = 0
            maxIter = 0
        end
    end

    for k = 1:maxIter
        type::String = "A="
        reset!(searchtracker)
        reset!(assemblytracker)

        v = t * G + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        accnorm = starnorm(v, I.dx)

        if accnorm <= S.beta
            if lastAccept <= S.kappa_updates[1]
                kappa = min(S.kappa, kappa^(S.kappa_powers[1]))
                type = "A+"
            elseif lastAccept >= S.kappa_updates[2]
                kappa = kappa^(S.kappa_powers[2])
                type = "A-"
            end
            lastAccept = 0
            lastAcceptedt = t
            lastAcceptedx = I.x

            t = min(t / kappa, t - S.gamma / starnorm(G, I, S.solveLS))
            I.dx = S.solveLS(I.hessianF, t * G + I.gradientF, I.P)
        else
            lastAccept += 1
            type = "S"
        end

        if lastAccept < S.kappa_updates[3]
            apply_descent!(I, S, tracker=searchtracker, assemblytracker=assemblytracker)
            handle_assembly!(A, L, assemblytracker,"A", k)
        else
            lastAccept = 0
            type = "R"
            t = lastAcceptedt
            kappa = kappa^(S.kappa_powers[3])
            set!(I, lastAcceptedx, S)
            
            if kappa < 1.000001
                handle_kappavanish!(A, "A", k)
                iterationCount = -k
                break
            end
        end

        critnorm = Inf
        if assemblytracker.singularity
            iterationCount = -k
            break
        else
            critnorm = starnorm(I.gradientF, I, S.solveLS)
        end

        log_iteration(
            L,
            k, 0 , "A",
            type, lastAccept, kappa, accnorm,
            searchtracker.i, searchtracker.val,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if critnorm <= bound
            iterationCount = k
            break
        end
        
        if k == S.maxIter
            handle_maxiterations!(A, "A", k)
            iterationCount = -k
            break
        end
    end

    if iterationCount > 0
        I.dx = S.solveLS(I.hessianF, I.gradientF, I.P)
        apply_descent!(I, S, useBacktracking=false, assemblytracker=assemblytracker)
        handle_assembly!(A, L, assemblytracker,"A", -1)

        if assemblytracker.singularity
            iterationCount *= -1
        else
            critnorm = starnorm(I.gradientF, I, S.solveLS)

            log_iteration(
                L,
                missing, 0 , "A",
                missing, missing, missing, missing,
                missing, missing,
                assemblytracker.conditionnumber,
                critnorm, missing, S.beta
            )

            if critnorm > S.beta
                handle_auxfail!(A)
                iterationCount *= -1
            end
        end
    end

    log_footer(L)

    A.Naux = iterationCount
end

"""
$(TYPEDSIGNATURES)

Executes main path-following with adaptive stepsize.

The iteration is performed on the [IterationData](@ref), which will later store the final
result.
If the iteration did not converge, the last iterate will still be provided as a result and 
the obtained accuracy stored in [AlgorithmData](@ref). 
The number of required iterations and potential messages will be stored in
[AlgorithmData](@ref) in any case.
If [LogData](@ref) is verbose, data per iteration will be written to the output stream.
Further, if a file is provided, the log will be also exported to that file.
"""
function pathfollowing_main_adaptive!(
    I::IterationData,
    A::AlgorithmData,
    S::StaticData,
    L::LogData
)
    searchtracker = DescentTracker(0,0.0)
    assemblytracker = AssemblyTracker(trackcondition=L.trackcondition)
    log_inital(L)

    kappa::Float64 = S.kappa
    t::Float64 = 0

    lastAccept::Int64 = 0
    lastAcceptedt::Float64 = t
    lastAcceptedx::AbstractVector{Float64} = I.x
    
    iterationCount::Int64 = 0
    
    for k = 1:S.maxIter
        type::String = "A="
        reset!(searchtracker)
        reset!(assemblytracker)

        v = t * S.c + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        accnorm = starnorm(v, I.dx)

        if accnorm <= S.beta
            if lastAccept <= S.kappa_updates[1]
                kappa = min(S.kappa, kappa^(S.kappa_powers[1]))
                type = "A+"
            elseif lastAccept >= S.kappa_updates[2]
                kappa = kappa^(S.kappa_powers[2])
                type = "A-"
            end
            lastAccept = 0
            lastAcceptedt = t
            lastAcceptedx = I.x

            t = max(kappa * t, t + (S.gamma / starnorm(S.c, I, S.solveLS)))
            I.dx = S.solveLS(I.hessianF, t * S.c + I.gradientF, I.P)
        else
            lastAccept += 1
            type = "S"
        end

        if lastAccept < S.kappa_updates[3]
            apply_descent!(I, S, tracker=searchtracker, assemblytracker=assemblytracker)
            handle_assembly!(A, L, assemblytracker,"M", k)

            if assemblytracker.singularity
                handle_accuracy!(A, S.tolFactor/t)
                iterationCount = -k
                break
            end
        else
            lastAccept = 0
            type = "R"
            t = lastAcceptedt
            kappa = kappa^(S.kappa_powers[3])
            set!(I, lastAcceptedx, S)
            
            if kappa < 1.000001
                handle_kappavanish!(A, "M", k)
                iterationCount = -k
                break
            end
        end

        log_iteration(
            L,
            k, A.Naux , "M",
            type, lastAccept, kappa, accnorm,
            searchtracker.i, searchtracker.val,
            assemblytracker.conditionnumber,
            missing, t, S.tolInv
        )
        
        if t >= S.tolInv
            iterationCount = k
            A.solution = I.x[1:S.lengthu]
            break
        end

        if k == S.maxIter
            handle_maxiterations!(A, "M", k)
            iterationCount = -k
            break
        end
    end

    log_footer(L)

    A.Nmain = iterationCount 
end

"""
$(TYPEDSIGNATURES)

Executes auxilliary path-following with long stepsize.

The iteration is performed on the [IterationData](@ref), which will later store
the final iterate as well as all the corresponding barrier terms.
The number of required iterations and potential messages will be stored in
[AlgorithmData](@ref).
If [LogData](@ref) is verbose, data per iteration will be written to the output stream.
Further, if a file is provided, the log will be also exported to that file.
"""
function pathfollowing_auxilliary_long!(
    I::IterationData,
    A::AlgorithmData,
    S::StaticData,
    L::LogData
)
    searchtracker = DescentTracker(0,0.0)
    assemblytracker = AssemblyTracker(trackcondition=L.trackcondition)
    log_inital(L)

    assemble!(I, S, tracker=assemblytracker)
    handle_assembly!(A, L, assemblytracker,"A", 0)
    
    t::Float64 = 1
    maxIter::Int64 = S.maxIter

    iterationCount::Int64 = 0
    G::AbstractVector{Float64} = -I.gradientF
    bound::Float64 = sqrt(S.beta) / (1 + sqrt(S.beta))

    if assemblytracker.singularity
        iterationCount = 0
        maxIter = 0
    else
        critnorm = starnorm(G, I, S.solveLS)

        log_iteration(
            L,
            0, 0 , "A",
            missing, missing, missing, missing,
            missing, missing,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if critnorm <= bound
            iterationCount = 0
            maxIter = 0
        end
    end

    for k = 1:maxIter
        type::String = "A"
        reset!(searchtracker)
        reset!(assemblytracker)
        
        v = t * G + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        accnorm = starnorm(v, I.dx)

        if accnorm <= S.beta
            t = min(t / S.kappa, t - S.gamma / starnorm(G, I, S.solveLS))
            I.dx = S.solveLS(I.hessianF, t * G + I.gradientF, I.P)
        else
            type = "S"
        end

        apply_descent!(I, S, tracker=searchtracker, assemblytracker=assemblytracker)
        handle_assembly!(A, L, assemblytracker,"A", k)  

        critnorm = Inf
        if assemblytracker.singularity
            iterationCount = -k
            break
        else
            critnorm = starnorm(I.gradientF, I, S.solveLS)
        end

        log_iteration(
            L,
            k, 0 , "A",
            type, missing, missing, accnorm,
            searchtracker.i, searchtracker.val,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if critnorm <= bound
            iterationCount = k
            break
        end
        
        if k == S.maxIter
            handle_maxiterations!(A, "A", k)
            iterationCount = -k
        end
    end

    if iterationCount > 0
        I.dx = S.solveLS(I.hessianF, I.gradientF, I.P)
        apply_descent!(I, S, useBacktracking=false, assemblytracker=assemblytracker)
        handle_assembly!(A, L, assemblytracker,"A", -1)

        if assemblytracker.singularity
            iterationCount *= -1
        else
            critnorm = starnorm(I.gradientF, I, S.solveLS)

            log_iteration(
                L,
                missing, 0 , "A",
                missing, missing, missing, missing,
                missing, missing,
                assemblytracker.conditionnumber,
                critnorm, missing, S.beta
            )

            if critnorm > S.beta
                handle_auxfail!(A)
                iterationCount *= -1
            end
        end
    end

    log_footer(L)

    A.Naux = iterationCount 
end

"""
$(TYPEDSIGNATURES)

Executes main path-following with long stepsize.

The iteration is performed on the [IterationData](@ref), which will later store the final
result.
If the iteration did not converge, the last iterate will still be provided as a result and 
the obtained accuracy stored in [AlgorithmData](@ref). 
The number of required iterations and potential messages will be stored in
[AlgorithmData](@ref) in any case.
If [LogData](@ref) is verbose, data per iteration will be written to the output stream.
Further, if a file is provided, the log will be also exported to that file.
"""
function pathfollowing_main_long!(
    I::IterationData,
    A::AlgorithmData,
    S::StaticData,
    L::LogData
)
    searchtracker = DescentTracker(0,0.0)
    assemblytracker = AssemblyTracker(trackcondition=L.trackcondition)
    log_inital(L)
    
    t::Float64 = 0
    iterationCount::Int64 = 0

    for k = 1:S.maxIter
        type::String = "A"
        reset!(searchtracker)
        reset!(assemblytracker)

        v = t * S.c + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        accnorm = starnorm(v, I.dx)


        if accnorm <= S.beta
            t = max(S.kappa * t, t + (S.gamma / starnorm(S.c, I, S.solveLS)))
            I.dx = S.solveLS(I.hessianF, t * S.c + I.gradientF, I.P)
        else
            type = "S"
        end

        apply_descent!(I, S, tracker=searchtracker, assemblytracker=assemblytracker)
        handle_assembly!(A, L, assemblytracker,"M", k)

        if assemblytracker.singularity
            handle_accuracy!(A, S.tolFactor/t)
            iterationCount = -k
            break
        end

        log_iteration(
            L,
            k, A.Naux , "M",
            type, missing, missing, accnorm,
            searchtracker.i, searchtracker.val,
            assemblytracker.conditionnumber,
            missing, t, S.tolInv
        )

        if t >= S.tolInv
            iterationCount = k
            A.solution = I.x[1:S.lengthu]
            break
        end

        if k == S.maxIter
            handle_maxiterations!(A, "M", k)
            iterationCount = -k
            break
        end
    end

    log_footer(L)

    A.Nmain = iterationCount 
end

"""
$(TYPEDSIGNATURES)

Executes auxilliary path-following with short stepsize.

The iteration is performed on the [IterationData](@ref), which will later store
the final iterate as well as all the corresponding barrier terms.
The number of required iterations and potential messages will be stored in
[AlgorithmData](@ref).
If [LogData](@ref) is verbose, data per iteration will be written to the output stream.
Further, if a file is provided, the log will be also exported to that file.
"""
function pathfollowing_auxilliary_short!(
    I::IterationData,
    A::AlgorithmData,
    S::StaticData,
    L::LogData
)
    assemblytracker = AssemblyTracker(trackcondition=L.trackcondition)
    log_inital(L)

    assemble!(I, S, tracker=assemblytracker)
    handle_assembly!(A, L, assemblytracker,"A", 0)
    
    t::Float64 = 1
    maxIter::Int64 = S.maxIter

    iterationCount::Int64 = 0
    G::AbstractVector{Float64} = -I.gradientF
    bound::Float64 = sqrt(S.beta) / (1 + sqrt(S.beta))

    if assemblytracker.singularity
        iterationCount = 0
        maxIter = 0
    else
        critnorm = starnorm(G, I, S.solveLS)

        log_iteration(
            L,
            0, 0 , "A",
            missing, missing, missing, missing,
            missing, missing,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if critnorm <= bound
            iterationCount = 0
            maxIter = 0
        end
    end
    
    for k = 1:maxIter
        reset!(assemblytracker)
        
        t -= S.gamma / starnorm(G, I, S.solveLS)
        I.dx = S.solveLS(I.hessianF, t * G + I.gradientF, I.P)
        apply_descent!(I, S, useBacktracking=false, assemblytracker=assemblytracker)
        handle_assembly!(A, L, assemblytracker,"A", k)

        critnorm = Inf
        if assemblytracker.singularity
            iterationCount = -k
            break
        else
            critnorm = starnorm(I.gradientF, I, S.solveLS)
        end
        
        log_iteration(
            L,
            k, 0 , "A",
            missing, missing, missing, missing,
            missing, missing,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )
        
        if critnorm <= bound
            iterationCount = k
            break
        end    
        
        if k == S.maxIter
            handle_maxiterations!(A, "A", k)
            iterationCount = -k
        end
    end
    
    if iterationCount > 0
        I.dx = S.solveLS(I.hessianF, I.gradientF, I.P)
        apply_descent!(I, S, useBacktracking=false, assemblytracker=assemblytracker)
        handle_assembly!(A, L, assemblytracker,"A", -1)

        if assemblytracker.singularity
            iterationCount *= -1
        else
            critnorm = starnorm(I.gradientF, I, S.solveLS)

            log_iteration(
                L,
                missing, 0 , "A",
                missing, missing, missing, missing,
                missing, missing,
                assemblytracker.conditionnumber,
                critnorm, missing, S.beta
            )

            if critnorm > S.beta
                handle_auxfail!(A)
                iterationCount *= -1
            end
        end
    end

    log_footer(L)

    A.Naux = iterationCount 
end

"""
$(TYPEDSIGNATURES)

Executes main path-following with short stepsize.

The iteration is performed on the [IterationData](@ref), which will later store the final
result.
If the iteration did not converge, the last iterate will still be provided as a result and 
the obtained accuracy stored in [AlgorithmData](@ref). 
The number of required iterations and potential messages will be stored in
[AlgorithmData](@ref) in any case.
If [LogData](@ref) is verbose, data per iteration will be written to the output stream.
Further, if a file is provided, the log will be also exported to that file.
"""
function pathfollowing_main_short!(
    I::IterationData,
    A::AlgorithmData,
    S::StaticData,
    L::LogData
)
    assemblytracker = AssemblyTracker(trackcondition=L.trackcondition)
    log_inital(L)

    t::Float64 = 0
    iterationCount::Int64 = 0
    
    for k = 1:S.maxIter
        reset!(assemblytracker)

        t += S.gamma / starnorm(S.c, I, S.solveLS)
        I.dx = S.solveLS(I.hessianF, t * S.c + I.gradientF, I.P)
        
        apply_descent!(I, S, useBacktracking=false, assemblytracker=assemblytracker)
        handle_assembly!(A, L, assemblytracker,"M", k)

        if assemblytracker.singularity
            handle_accuracy!(A, S.tolFactor/t)
            iterationCount = -k
            break
        end

        log_iteration(
            L,
            k, A.Naux , "M",
            missing, missing, missing, missing,
            missing, missing,
            assemblytracker.conditionnumber,
            missing, t, S.tolInv
        )
        
        if t >= S.tolInv
            iterationCount = k
            A.solution = I.x[1:S.lengthu]
            break
        end

        if k == S.maxIter
            handle_maxiterations!(A, "M", k)
            iterationCount = -k
            break
        end
    end

    log_footer(L)

    A.Nmain = iterationCount 
end

"""
    select_pathfollowing(stepsize::Stepsize) -> Tuple{Function, Function}
    
Returns functions for auxilliary and main pathfollowing dependend on the stepsize.
"""
function select_pathfollowing(stepsize::Stepsize) :: Tuple{Function, Function}
    if stepsize === LONG
        return pathfollowing_auxilliary_long!, pathfollowing_main_long!
    elseif stepsize === ADAPTIVE
        return pathfollowing_auxilliary_adaptive!, pathfollowing_main_adaptive!
    else
        return pathfollowing_auxilliary_short!, pathfollowing_main_short!
    end
end
