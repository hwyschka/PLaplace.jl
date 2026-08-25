"""
$(TYPEDEF)

Public type used to specify the stepping scheme of the path-following schemes.
For more information see the [interior-point section](@ref path-following-theory)
in the documentation.

# Available Options
- `SHORT`:
    Regular update of the parameter t.
- `LONG`:
    Larger update of the parameter t if iterate fulfill approximate centering condition.
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

Base.parse(E::Type{<:Enum}, str::String) =
    let insts = instances(E) ,
        p = findfirst(==(Symbol(str)) ∘ Symbol, insts) ;
        p !== nothing ? insts[p] : nothing
    end

"""
$(TYPEDSIGNATURES)

Executes auxiliary path-following with adaptive stepsize.

The iteration is performed on the [IterationData](@ref), which will later store
the final iterate as well as all the corresponding barrier terms.
The number of required iterations and potential messages will be stored in
[AlgorithmData](@ref).
If [LogData](@ref) is verbose, data per iteration will be written to the output stream.
Further, if a file is provided, the log will be also exported to that file.
"""
function pathfollowing_auxiliary_adaptive!(
    I::IterationData,
    A::AlgorithmData,
    S::StaticData,
    L::LogData
)
    starnormtracker = StarnormTracker()
    descenttracker = DescentTracker()
    assemblytracker = AssemblyTracker(
        trackcondition = L.trackcondition,
        trackhessian = L.exporthessian
    )
    log_inital(L)

    assemble!(I, S, tracker=assemblytracker)
    handle_assembly!(A, L, assemblytracker, "A", 0)
    
    t::Float64 = 1
    kappa::Float64 = S.kappa
    maxIter::Int64 = S.maxIter
    
    iterationCount::Int64 = 0
    G::AbstractVector{Float64} = -I.gradientF
    bound::Float64 = sqrt(S.beta) / (1 + sqrt(S.beta))

    lastAccept::Int64 = 0
    lastAcceptedt::Float64 = t
    lastAcceptedx::AbstractVector{Float64} = I.x

    rejection::Bool = false
    isshort::Bool = false
    skipfinalupdate::Bool = false

    if assemblytracker.singularity
        iterationCount = -1
        maxIter = 0
    else
        critnorm = starnorm(G, I, S.solveLS, tracker=starnormtracker)

        log_iteration(
            L,
            0, 0 , "A",
            missing, missing, missing, missing,
            missing, missing,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if starnormtracker.failed
            handle_starnorm!(A, "crit", "A", 0)
            iterationCount = -1
            maxIter = 0
            skipfinalupdate = true
        end

        if critnorm <= S.beta
            iterationCount = 0
            maxIter = 0
            skipfinalupdate = true
        elseif critnorm <= bound
            iterationCount = 0
            maxIter = 0
        end
    end

    for k = 1:maxIter
        type::String = "A="
        reset!(descenttracker)
        reset!(assemblytracker)

        v = t * G + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        accnorm = starnorm(v, I.dx, tracker=starnormtracker)

        if starnormtracker.failed
            handle_starnorm!(A, "acc", "A", k)
            iterationCount = -k
            break
        end

        if accnorm <= S.beta
            if rejection
                if isshort
                    handle_kappavanish!(A, "A", k)
                    iterationCount = -k
                    break
                end

                kappa = kappa^(S.kappa_powers[3])
                type = "R"
            elseif lastAccept >= S.kappa_updates[2]
                kappa = kappa^(S.kappa_powers[2])
                type = "A-"
            elseif lastAccept <= S.kappa_updates[1]
                kappa = min(S.kappa, kappa^(S.kappa_powers[1]))
                type = "A+"
            end
            lastAccept = 0
            lastAcceptedt = t
            lastAcceptedx = I.x

            Gnorm = starnorm(G, I, S.solveLS, tracker=starnormtracker)

            if starnormtracker.failed
                handle_starnorm!(A, "G", "A", k)
                iterationCount = -k
                break
            end

            t_short = t - S.gamma / Gnorm
            t_long = t / kappa
            if t_long <  t_short
                t = t_long
                isshort = false
            else
                t = t_short
                isshort = true
            end

            v = t * G + I.gradientF
            I.dx = S.solveLS(I.hessianF, v, I.P)
        else
            lastAccept += 1
            type = "S"
        end

        if lastAccept < S.kappa_updates[3]
            rejection = false
            apply_descent!(
                I, t, G, v, S,
                tracker=descenttracker, assemblytracker=assemblytracker
            )
            handle_descent!(A, L, descenttracker, assemblytracker, "A", k)
        else
            rejection = true
            t = lastAcceptedt
            set!(I, lastAcceptedx, S)
        end

        if descenttracker.failed || assemblytracker.singularity
            log_iteration(
                L,
                k, 0 , "A",
                type, lastAccept, kappa, accnorm,
                descenttracker.i, descenttracker.val,
                assemblytracker.conditionnumber,
                missing, t, bound
            )

            iterationCount = -k
            break
        end
        
        critnorm = starnorm(I.gradientF, I, S.solveLS, tracker=starnormtracker)

        log_iteration(
            L,
            k, 0 , "A",
            type, lastAccept, kappa, accnorm,
            descenttracker.i, descenttracker.val,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if starnormtracker.failed
            handle_starnorm!(A, "crit", "A", k)
            iterationCount = -k
            break
        end

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

    if iterationCount >= 0 && !skipfinalupdate
        I.dx = S.solveLS(I.hessianF, I.gradientF, I.P)
        apply_descent!(
            I, missing, missing, I.gradientF, S,
            force_nobacktracking=true,
            tracker=descenttracker, assemblytracker=assemblytracker
        )
        handle_descent!(A, L, descenttracker, assemblytracker, "A", -1)

        if descenttracker.failed || assemblytracker.singularity
            iterationCount = iterationCount > 0 ? -iterationCount : -1
        else
            critnorm = starnorm(I.gradientF, I, S.solveLS, tracker=starnormtracker)

            log_iteration(
                L,
                missing, 0 , "A",
                missing, missing, missing, missing,
                missing, missing,
                assemblytracker.conditionnumber,
                critnorm, missing, S.beta
            )

            if starnormtracker.failed
                handle_starnorm!(A, "crit", "A", -1)
                iterationCount = iterationCount > 0 ? -iterationCount : -1
            end

            if critnorm > S.beta
                handle_auxfail!(A)
                iterationCount = iterationCount > 0 ? -iterationCount : -1
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
    starnormtracker = StarnormTracker()
    descenttracker = DescentTracker()
    assemblytracker = AssemblyTracker(
        trackcondition = L.trackcondition,
        trackhessian = L.exporthessian
    )
    log_inital(L)

    kappa::Float64 = S.kappa
    t::Float64 = 0

    lastAccept::Int64 = 0
    lastAcceptedt::Float64 = t
    lastAcceptedx::AbstractVector{Float64} = I.x

    rejection::Bool = false # leads to automatic acceptance of first step by construction
    isshort::Bool = false
    
    iterationCount::Int64 = 0
    
    for k = 1:S.maxIter
        type::String = "A="
        reset!(descenttracker)
        reset!(assemblytracker)
 
        v = t * S.c + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        accnorm = starnorm(v, I.dx, tracker=starnormtracker)

        if starnormtracker.failed
            handle_starnorm!(A, "acc", "M", k)
            iterationCount = -k
            break
        end

        if accnorm <= S.beta
            if t >= S.tolInv
                iterationCount = k
                A.solution = I.x[1:S.lengthu]
                break
            end

            if rejection
                if isshort
                    handle_kappavanish!(A, "M", k)
                    iterationCount = -k
                    break
                end

                kappa = kappa^(S.kappa_powers[3])
                type = "R"
                rejection = false
            elseif lastAccept >= S.kappa_updates[2]
                kappa = kappa^(S.kappa_powers[2])
                type = "A-"
            elseif lastAccept <= S.kappa_updates[1]
                kappa = min(S.kappa, kappa^(S.kappa_powers[1]))
                type = "A+"
            end
            lastAccept = 0
            lastAcceptedt = t
            lastAcceptedx = I.x

            cnorm = starnorm(S.c, I, S.solveLS, tracker=starnormtracker)

            if starnormtracker.failed
                handle_starnorm!(A, "c", "M", k)
                iterationCount = -k
                break
            end

            t_short = t + (S.gamma / cnorm)
            t_long = t * kappa
            if t_long >  t_short
                t = t_long
                isshort = false
            else
                t = t_short
                isshort = true
            end

            v = t * S.c + I.gradientF
            I.dx = S.solveLS(I.hessianF, v, I.P)
        else
            lastAccept += 1
            type = "S"
        end

        if lastAccept < S.kappa_updates[3]
            apply_descent!(
                I, t, S.c, v, S,
                tracker=descenttracker, assemblytracker=assemblytracker
            )
            handle_descent!(A, L, descenttracker, assemblytracker, "M", k)

            log_iteration(
                L,
                k, A.Naux , "M",
                type, lastAccept, kappa, accnorm,
                descenttracker.i, descenttracker.val,
                assemblytracker.conditionnumber,
                missing, t, S.tolInv
            )

            if descenttracker.failed || assemblytracker.singularity
                handle_accuracy!(A, S.tolFactor/t)
                iterationCount = -k
                break
            end
        else
            log_iteration(
                L,
                k, A.Naux , "M",
                type, lastAccept, kappa, accnorm,
                missing, missing,
                missing,
                missing, t, S.tolInv
            )
            
            rejection = true
            t = lastAcceptedt
            set!(I, lastAcceptedx, S)
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

Executes auxiliary path-following with long stepsize.

The iteration is performed on the [IterationData](@ref), which will later store
the final iterate as well as all the corresponding barrier terms.
The number of required iterations and potential messages will be stored in
[AlgorithmData](@ref).
If [LogData](@ref) is verbose, data per iteration will be written to the output stream.
Further, if a file is provided, the log will be also exported to that file.
"""
function pathfollowing_auxiliary_long!(
    I::IterationData,
    A::AlgorithmData,
    S::StaticData,
    L::LogData
)
    starnormtracker = StarnormTracker()
    descenttracker = DescentTracker()
    assemblytracker = AssemblyTracker(
        trackcondition = L.trackcondition,
        trackhessian = L.exporthessian
    )
    log_inital(L)

    assemble!(I, S, tracker=assemblytracker)
    handle_assembly!(A, L, assemblytracker, "A", 0)
    
    t::Float64 = 1
    maxIter::Int64 = S.maxIter

    iterationCount::Int64 = 0
    G::AbstractVector{Float64} = -I.gradientF
    bound::Float64 = sqrt(S.beta) / (1 + sqrt(S.beta))

    skipfinalupdate::Bool = false

    if assemblytracker.singularity
        iterationCount = -1
        maxIter = 0
    else
        critnorm = starnorm(G, I, S.solveLS, tracker=starnormtracker)

        log_iteration(
            L,
            0, 0 , "A",
            missing, missing, missing, missing,
            missing, missing,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if starnormtracker.failed
            handle_starnorm!(A, "crit", "A", 0)
            iterationCount = -1
            maxIter = 0
            skipfinalupdate = true
        end

        if critnorm <= S.beta
            iterationCount = 0
            maxIter = 0
            skipfinalupdate = true
        elseif critnorm <= bound
            iterationCount = 0
            maxIter = 0
        end
    end

    for k = 1:maxIter
        type::String = "A"
        reset!(descenttracker)
        reset!(assemblytracker)
        
        v = t * G + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        accnorm = starnorm(v, I.dx, tracker=starnormtracker)

        if starnormtracker.failed
            handle_starnorm!(A, "acc", "A", k)
            iterationCount = -k
            break
        end

        if accnorm <= S.beta
            Gnorm = starnorm(G, I, S.solveLS, tracker=starnormtracker)

            if starnormtracker.failed
                handle_starnorm!(A, "G", "A", k)
                iterationCount = -k
                break
            end

            t = min(t / S.kappa, t - S.gamma / Gnorm)
            v = t * G + I.gradientF
            I.dx = S.solveLS(I.hessianF, v, I.P)
        else
            type = "S"
        end

        apply_descent!(
            I, t, G, v, S,
            tracker=descenttracker, assemblytracker=assemblytracker
        )
        handle_descent!(A, L, descenttracker, assemblytracker, "A", k)  

        if descenttracker.failed || assemblytracker.singularity
            log_iteration(
                L,
                k, 0 , "A",
                type, lastAccept, kappa, accnorm,
                descenttracker.i, descenttracker.val,
                assemblytracker.conditionnumber,
                missing, t, bound
            )

            iterationCount = -k
            break
        end
        
        critnorm = starnorm(I.gradientF, I, S.solveLS, tracker=starnormtracker)

        log_iteration(
            L,
            k, 0 , "A",
            type, missing, missing, accnorm,
            descenttracker.i, descenttracker.val,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if starnormtracker.failed
            handle_starnorm!(A, "crit", "A", k)
            iterationCount = -k
            break
        end

        if critnorm <= bound
            iterationCount = k
            break
        end
        
        if k == S.maxIter
            handle_maxiterations!(A, "A", k)
            iterationCount = -k
        end
    end

    if iterationCount >= 0 && !skipfinalupdate
        I.dx = S.solveLS(I.hessianF, I.gradientF, I.P)
        apply_descent!(
            I, missing, missing, I.gradientF, S,
            force_nobacktracking=true,
            tracker=descenttracker, assemblytracker=assemblytracker
        )
        handle_descent!(A, L, descenttracker, assemblytracker, "A", -1)

        if descenttracker.failed || assemblytracker.singularity
            iterationCount = iterationCount > 0 ? -iterationCount : -1
        else
            critnorm = starnorm(I.gradientF, I, S.solveLS, tracker=starnormtracker)

            log_iteration(
                L,
                missing, 0 , "A",
                missing, missing, missing, missing,
                missing, missing,
                assemblytracker.conditionnumber,
                critnorm, missing, S.beta
            )

            if starnormtracker.failed
                handle_starnorm!(A, "crit", "A", -1)
                iterationCount = iterationCount > 0 ? -iterationCount : -1
            end

            if critnorm > S.beta
                handle_auxfail!(A)
                iterationCount = iterationCount > 0 ? -iterationCount : -1
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
    starnormtracker = StarnormTracker()
    descenttracker = DescentTracker()
    assemblytracker = AssemblyTracker(
        trackcondition = L.trackcondition,
        trackhessian = L.exporthessian
    )
    log_inital(L)
    
    t::Float64 = 0
    iterationCount::Int64 = 0

    for k = 1:S.maxIter
        type::String = "A"
        reset!(descenttracker)
        reset!(assemblytracker)

        v = t * S.c + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        accnorm = starnorm(v, I.dx, tracker=starnormtracker)

        if starnormtracker.failed
            handle_starnorm!(A, "acc", "M", k)
            iterationCount = -k
            break
        end

        if accnorm <= S.beta
            if t >= S.tolInv
                iterationCount = k
                A.solution = I.x[1:S.lengthu]
                break
            end

            cnorm = starnorm(S.c, I, S.solveLS, tracker=starnormtracker)

            if starnormtracker.failed
                handle_starnorm!(A, "c", "M", k)
                iterationCount = -k
                break
            end

            t = max(S.kappa * t, t + (S.gamma / cnorm))
            v = t * S.c + I.gradientF
            I.dx = S.solveLS(I.hessianF, v, I.P)
        else
            type = "S"
        end

        apply_descent!(
            I, t, S.c, v, S,
            tracker=descenttracker, assemblytracker=assemblytracker
        )
        handle_descent!(A, L, descenttracker, assemblytracker, "M", k)

        log_iteration(
            L,
            k, A.Naux , "M",
            type, missing, missing, accnorm,
            descenttracker.i, descenttracker.val,
            assemblytracker.conditionnumber,
            missing, t, S.tolInv
        )
        
        if descenttracker.failed || assemblytracker.singularity
            handle_accuracy!(A, S.tolFactor/t)
            iterationCount = -k
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

Executes auxiliary path-following with short stepsize.

The iteration is performed on the [IterationData](@ref), which will later store
the final iterate as well as all the corresponding barrier terms.
The number of required iterations and potential messages will be stored in
[AlgorithmData](@ref).
If [LogData](@ref) is verbose, data per iteration will be written to the output stream.
Further, if a file is provided, the log will be also exported to that file.
"""
function pathfollowing_auxiliary_short!(
    I::IterationData,
    A::AlgorithmData,
    S::StaticData,
    L::LogData
)
    starnormtracker = StarnormTracker()
    descenttracker = DescentTracker()
    assemblytracker = AssemblyTracker(
        trackcondition = L.trackcondition,
        trackhessian = L.exporthessian
    )
    log_inital(L)

    assemble!(I, S, tracker=assemblytracker)
    handle_assembly!(A, L, assemblytracker,"A", 0)
    
    t::Float64 = 1
    maxIter::Int64 = S.maxIter

    iterationCount::Int64 = 0
    G::AbstractVector{Float64} = -I.gradientF
    bound::Float64 = sqrt(S.beta) / (1 + sqrt(S.beta))

    skipfinalupdate::Bool = false

    if assemblytracker.singularity
        iterationCount = -1
        maxIter = 0
    else
        critnorm = starnorm(G, I, S.solveLS, tracker=starnormtracker)

        log_iteration(
            L,
            0, 0 , "A",
            missing, missing, missing, missing,
            missing, missing,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if starnormtracker.failed
            handle_starnorm!(A, "crit", "A", 0)
            iterationCount = -1
            maxIter = 0
            skipfinalupdate = true
        end

        if critnorm <= S.beta
            iterationCount = 0
            maxIter = 0
            skipfinalupdate = true
        elseif critnorm <= bound
            iterationCount = 0
            maxIter = 0
        end
    end
    
    for k = 1:maxIter
        reset!(descenttracker)
        reset!(assemblytracker)
        
        Gnorm = starnorm(G, I, S.solveLS, tracker=starnormtracker)

        if starnormtracker.failed
            handle_starnorm!(A, "G", "A", k)
            iterationCount = -k
            break
        end

        t -= S.gamma / Gnorm
        v = t * G + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        apply_descent!(
            I, t, G, v, S,
            backtracking=false, 
            tracker=descenttracker, assemblytracker=assemblytracker
        )
        handle_descent!(A, L, descenttracker, assemblytracker, "A", k)

        if descenttracker.failed || assemblytracker.singularity
            log_iteration(
                L,
                k, 0 , "A",
                missing, missing, missing, missing,
                descenttracker.i, descenttracker.val,
                assemblytracker.conditionnumber,
                missing, t, bound
            )

            iterationCount = -k
            break
        end
            
        critnorm = starnorm(I.gradientF, I, S.solveLS, tracker=starnormtracker)

        log_iteration(
            L,
            k, 0 , "A",
            missing, missing, missing, missing,
            descenttracker.i, descenttracker.val,
            assemblytracker.conditionnumber,
            critnorm, t, bound
        )

        if starnormtracker.failed
            handle_starnorm!(A, "crit", "A", k)
            iterationCount = -k
            break
        end
        
        if critnorm <= bound
            iterationCount = k
            break
        end    
        
        if k == S.maxIter
            handle_maxiterations!(A, "A", k)
            iterationCount = -k
        end
    end
    
    if iterationCount >= 0 && !skipfinalupdate
        I.dx = S.solveLS(I.hessianF, I.gradientF, I.P)
        apply_descent!(
            I, missing, missing, I.gradientF, S,
            force_nobacktracking=true,
            tracker=descenttracker, assemblytracker=assemblytracker 
        )
        handle_descent!(A, L, descenttracker, assemblytracker, "A", -1)

        if descenttracker.failed || assemblytracker.singularity
            iterationCount = iterationCount > 0 ? -iterationCount : -1
        else
            critnorm = starnorm(I.gradientF, I, S.solveLS, tracker=starnormtracker)

            log_iteration(
                L,
                missing, 0 , "A",
                missing, missing, missing, missing,
                missing, missing,
                assemblytracker.conditionnumber,
                critnorm, missing, S.beta
            )

            if starnormtracker.failed
                handle_starnorm!(A, "crit", "A", -1)
                iterationCount = iterationCount > 0 ? -iterationCount : -1
            end

            if critnorm > S.beta
                handle_auxfail!(A)
                iterationCount = iterationCount > 0 ? -iterationCount : -1
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
    starnormtracker = StarnormTracker()
    descenttracker = DescentTracker()
    assemblytracker = AssemblyTracker(
        trackcondition = L.trackcondition,
        trackhessian = L.exporthessian
    )
    log_inital(L)

    t::Float64 = 0
    iterationCount::Int64 = 0
    
    for k = 1:S.maxIter
        reset!(descenttracker)
        reset!(assemblytracker)

        cnorm = starnorm(S.c, I, S.solveLS, tracker=starnormtracker)

        if starnormtracker.failed
            handle_starnorm!(A, "c", "M", k)
            iterationCount = -k
            break
        end

        t += S.gamma / cnorm
        v = t * S.c + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        
        apply_descent!(
            I, t, S.c, v, S,
            backtracking=false,
            tracker=descenttracker, assemblytracker=assemblytracker
        )
        handle_descent!(A, L, descenttracker, assemblytracker, "M", k)

        log_iteration(
            L,
            k, A.Naux , "M",
            missing, missing, missing, missing,
            descenttracker.i, descenttracker.val,
            assemblytracker.conditionnumber,
            missing, t, S.tolInv
        )

        if descenttracker.failed || assemblytracker.singularity
            handle_accuracy!(A, S.tolFactor/t)
            iterationCount = -k
            break
        end
        
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
    
Returns functions for auxiliary and main pathfollowing dependend on the stepsize.
"""
function select_pathfollowing(stepsize::Stepsize) :: Tuple{Function, Function}
    if stepsize === LONG
        return pathfollowing_auxiliary_long!, pathfollowing_main_long!
    elseif stepsize === ADAPTIVE
        return pathfollowing_auxiliary_adaptive!, pathfollowing_main_adaptive!
    else
        return pathfollowing_auxiliary_short!, pathfollowing_main_short!
    end
end
