"""
    pathfollowing_auxilliary_adaptive(S::StaticData, I::IterationData)

Executes auxilliary path-following with adaptive stepsize on the iteration data I and returns required iterations.
"""
function pathfollowing_auxilliary_adaptive(S::StaticData, I::IterationData)
    t::Float64 = 1
    kappa::Float64 = S.kappa
    
    iterationCount::Int64 = 0
    G::AbstractVector{Float64} = -I.gradientF
    bound::Float64 = sqrt(S.beta)/(1+sqrt(S.beta))

    lastAccept::Int64 = 0
    lastAcceptedt::Float64 = t
    lastAcceptedx::AbstractVector{Float64} = I.x

    if S.logObjective
        write_step(S.auxObjectiveFileName, 0, starnorm(G, I, S.solveLS))
    end

    if S.logCondition
        write_step(S.auxConditionFileName, 0, cond(I.unfactorizedHessianF,Inf))
    end

    for k = 1:S.maxIter
        S.log("========================================================")
        S.log("Auxiliary path following scheme - Step $k")
        S.log("========================================================")
        
        v = t * G + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        snorm = starnorm(v, I.dx)

        if snorm <= S.beta 
            S.log("Acceptance: $snorm ≤ $(S.beta)")
            if lastAccept <= 2
                kappa = min(S.kappa, kappa^2)
            elseif lastAccept >= 8
                kappa = sqrt(kappa)
            end
            lastAccept = 0
            lastAcceptedt = t
            lastAcceptedx = I.x

            t = min(t / kappa, t - S.gamma / starnorm(G, I, S.solveLS))
            I.dx = S.solveLS(I.hessianF, t * G + I.gradientF, I.P)
        else
            S.log("slow step: $snorm > $(S.beta)")
            lastAccept += 1
        end

        if lastAccept < 15
            apply_descent!(I,S)
        else
            lastAccept = 0
            t = lastAcceptedt
            kappa = kappa^(1/4)
            set!(I, lastAcceptedx, S)
            
            if kappa < 1.000001
                iterationCount = -k
                I.msg *= " ⲕ too small  - "
                S.log(I.msg)
                break
            end
        end

        if S.logCondition
            write_step(S.auxConditionFileName, k, cond(I.unfactorizedHessianF,Inf))
        end

        try
            snorm = starnorm(I.gradientF, I, S.solveLS)
        catch e
            if isa(e, SingularException)
                S.log("Hessian singular at t = $t")
                I.msg *= " Hessian singular - "
                iterationCount = -k    
                break
            else
                rethrow()
            end
        end
        
        if S.logObjective
            write_step(S.auxObjectiveFileName, k, snorm)
        end

        if snorm <= bound
            S.log("Norm: $snorm ≤ $bound, t = $t")
            iterationCount = k
            break
        end
        S.log("Norm $snorm > $bound, t = $t")
        
        if k == S.maxIter
            S.log("Maximum number of iterations reached.")
            I.msg *= " Exceeded iterations  - "
            iterationCount = -k
        end
    end

    if iterationCount > 0
        I.dx = S.solveLS(I.hessianF, I.gradientF, I.P)
        apply_descent!(I, S, useBacktracking=false)
    end

    return iterationCount 
end

"""
    pathfollowing_main_adaptive(S::StaticData, I::IterationData)

Executes main path-following with adaptive stepsize on the iteration data I and returns required iterations.
"""
function pathfollowing_main_adaptive(S::StaticData, I::IterationData)
    kappa::Float64 = S.kappa
    t::Float64 = 0

    lastAccept::Int64 = 0
    lastAcceptedt::Float64 = t
    lastAcceptedx::AbstractVector{Float64} = I.x
    
    iterationCount::Int64 = 0

    if S.logCondition
        write_step(S.mainConditionFileName, 0, cond(I.unfactorizedHessianF,Inf))
    end
    
    for k = 1:S.maxIter
        S.log("========================================================")
        S.log("Main path following scheme - Step $k")
        S.log("========================================================")

        snorm::Float64 = 0
        try
            v = t * S.c + I.gradientF
            I.dx = S.solveLS(I.hessianF, v, I.P)
            snorm = starnorm(v, I.dx)
        catch e
            if isa(e, SingularException)
                S.log("Hessian singular at t = $t < $(S.tolInv)")
                I.msg *= " Hessian singular  - "
                iterationCount = -(k-1)
                S.eps = S.tolFactor / t    
                break
            else
                rethrow()
            end
        end

        if snorm <= S.beta
            S.log("t,x accepted: $snorm ≤ $(S.beta)")
            if lastAccept <= 2
                kappa = min(S.kappa, kappa^2)
            elseif lastAccept >= 8
                kappa = sqrt(kappa)
            end
            lastAccept = 0
            lastAcceptedt = t
            lastAcceptedx = I.x

            t = max(kappa * t, t + (S.gamma / starnorm(S.c, I, S.solveLS)))
            I.dx = S.solveLS(I.hessianF, t * S.c + I.gradientF, I.P)
        else
            S.log("slow step: $snorm > $(S.beta)")
            lastAccept += 1
        end

        if lastAccept < 15
            apply_descent!(I,S)
        else
            lastAccept = 0
            t = lastAcceptedt
            kappa = kappa^(1/4)
            set!(I, lastAcceptedx, S)
            
            S.log("t,x rejected. New ⲕ = $(kappa)")
            if kappa < 1.000001
                iterationCount = -k
                S.eps = S.tolFactor / t   
                I.msg *= " ⲕ too small  - "
                S.log(I.msg)
                break
            end
        end

        if S.logObjective
            write_step(S.mainObjectiveFileName, k, S.tolFactor / t)
        end

        if S.logCondition
            write_step(S.mainConditionFileName, k, cond(I.unfactorizedHessianF,Inf))
        end
        
        if t >= S.tolInv
            S.log("t = $t ≥ $(S.tolInv)")
            iterationCount = k
            break
        end
        S.log("t = $t < $(S.tolInv)")

        if k == S.maxIter
            S.log("Maximum number of iterations reached.")
            I.msg *= " Exceeded iterations  - "
            iterationCount = -k
            break
        end
    end

    return iterationCount
end

"""
    pathfollowing_auxilliary_long(S::StaticData, I::IterationData)

Executes auxilliary path-following with long stepsize on the iteration data I and returns required iterations.
"""
function pathfollowing_auxilliary_long(S::StaticData, I::IterationData)
    t::Float64 = 1

    iterationCount::Int64 = 0
    G::AbstractVector{Float64} = -I.gradientF
    bound::Float64 = sqrt(S.beta)/(1+sqrt(S.beta))

    if S.logObjective
        write_step(S.auxObjectiveFileName, 0, starnorm(G, I, S.solveLS))
    end

    if S.logCondition
        write_step(S.auxConditionFileName, 0, cond(I.unfactorizedHessianF,Inf))
    end

    for k = 1:S.maxIter
        S.log("========================================================")
        S.log("Auxiliary path following scheme - Step $k")
        S.log("========================================================")
        
        v = t * G + I.gradientF
        I.dx = S.solveLS(I.hessianF, v, I.P)
        snorm = starnorm(v, I.dx)

        if snorm <= S.beta
            S.log("t,x accepted: $snorm ≤ $(S.beta)") 
            t = min(t / S.kappa, t - S.gamma / starnorm(G, I, S.solveLS))
            I.dx = S.solveLS(I.hessianF, t * G + I.gradientF, I.P)
        else
            S.log("slow step: $snorm > $(S.beta)")
        end

        apply_descent!(I, S)    
        
        if S.logCondition
            write_step(S.auxConditionFileName, k, cond(I.unfactorizedHessianF,Inf))
        end

        try
            snorm = starnorm(I.gradientF, I, S.solveLS)
        catch e
            if isa(e, SingularException)
                S.log("Hessian singular at t = $t")
                I.msg *= " Hessian singular - "
                iterationCount = -k    
                break
            else
                rethrow()
            end
        end
        
        if S.logObjective
            write_step(S.auxObjectiveFileName, k, snorm)
        end

        if snorm <= bound
            S.log("Norm: $snorm ≤ $bound, t = $t")
            iterationCount = k
            break
        end
        S.log("Norm: $snorm > $bound, t = $t")
        
        if k == S.maxIter
            S.log("Maximum number of iterations reached.")
            I.msg *= " Exceeded iterations  - "
            iterationCount = -k
        end
    end

    if iterationCount > 0
        I.dx = S.solveLS(I.hessianF, I.gradientF, I.P)
        apply_descent!(I, S, useBacktracking=false)
    end

    return iterationCount 
end

"""
    pathfollowing_main_long(S::StaticData, I::IterationData)

Executes main path-following with long stepsize on the iteration data I and returns required iterations.
"""
function pathfollowing_main_long(S::StaticData, I::IterationData)
    t::Float64 = 0
    iterationCount::Int64 = 0

    if S.logCondition
        write_step(S.mainConditionFileName, 0, cond(I.unfactorizedHessianF,Inf))
    end

    for k = 1:S.maxIter
        S.log("========================================================")
        S.log("Main path following scheme - Step $k")
        S.log("========================================================")

        snorm::Float64 = 0
        try
            v = t * S.c + I.gradientF
            I.dx = S.solveLS(I.hessianF, v, I.P)
            snorm = starnorm(v, I.dx)
        catch e
            if isa(e, SingularException)
                S.log("Hessian singular at t = $t < $(S.tolInv)")
                I.msg *= " Hessian singular  - "
                iterationCount = -(k-1)
                S.eps = S.tolFactor / t    
                break
            else
                rethrow()
            end
        end

        if snorm <= S.beta
            S.log("t,x accepted: $snorm ≤ $(S.beta)") 
            t = max(S.kappa * t, t + (S.gamma / starnorm(S.c, I, S.solveLS)))
            I.dx = S.solveLS(I.hessianF, t * S.c + I.gradientF, I.P)
        else
            S.log("slow step: $snorm > $(S.beta)")
        end

        apply_descent!(I, S)

        if S.logObjective
            write_step(S.mainObjectiveFileName, k, S.tolFactor / t)
        end

        if S.logCondition
            write_step(S.mainConditionFileName, k, cond(I.unfactorizedHessianF,Inf))
        end

        if t >= S.tolInv
            S.log("t = $t ≥ $(S.tolInv)")
            iterationCount = k
            break
        end
        S.log("t = $t < $(S.tolInv)")

        if k == S.maxIter
            S.log("Maximum number of iterations reached.")
            I.msg *= " Exceeded iterations  - "
            iterationCount = -k
            break
        end
    end

    return iterationCount
end

"""
    pathfollowing_auxilliary_short(S::StaticData, I::IterationData)

Executes auxilliary path-following with short stepsize on the iteration data I and returns required iterations.
"""
function pathfollowing_auxilliary_short(S::StaticData, I::IterationData)
    t::Float64 = 1

    iterationCount::Int64 = 0
    G::AbstractVector{Float64} = -I.gradientF
    bound::Float64 = sqrt(S.beta)/(1+sqrt(S.beta))

    if S.logObjective
        write_step(S.auxObjectiveFileName, 0, starnorm(G, I, S.solveLS))
    end

    if S.logCondition
        write_step(S.auxConditionFileName, 0, cond(I.unfactorizedHessianF,Inf))
    end
    
    for k = 1:S.maxIter
        S.log("========================================================")
        S.log("Auxiliary path following scheme - Step $k")
        S.log("========================================================")
        
        t -= S.gamma / starnorm(G, I, S.solveLS)
        I.dx = S.solveLS(I.hessianF, t * G + I.gradientF, I.P)
        apply_descent!(I, S, useBacktracking=false)

        if S.logCondition
            write_step(S.auxConditionFileName, k, cond(I.unfactorizedHessianF,Inf))
        end

        snorm::Float64 = 0
        try
            snorm = starnorm(I.gradientF, I, S.solveLS)
        catch e
            if isa(e, SingularException)
                S.log("Hessian singular at t = $t")
                I.msg *= " Hessian singular - "
                iterationCount = -k    
                break
            else
                rethrow()
            end
        end
        
        if S.logObjective
            write_step(S.auxObjectiveFileName, k, snorm)
        end
        
        if snorm <= bound
            S.log("Norm: $snorm ≤ $bound, t = $t")
            iterationCount = k
            break
        end
        S.log("Norm: $snorm > $bound, t = $t")      
        
        if k == S.maxIter
            S.log("Maximum number of iterations reached.")
            I.msg *= " Exceeded iterations  - "
            iterationCount = -k
        end
    end
    
    if iterationCount > 0
        I.dx = S.solveLS(I.hessianF, I.gradientF, I.P)
        apply_descent!(I, S, useBacktracking=false)
    end

    return iterationCount 
end

"""
    pathfollowing_main_short(S::StaticData, I::IterationData)

Executes main path-following with short stepsize on the iteration data I and returns required iterations.
"""
function pathfollowing_main_short(S::StaticData, I::IterationData)
    t::Float64 = 0
    iterationCount::Int64 = 0

    if S.logCondition
        write_step(S.mainConditionFileName, 0, cond(I.unfactorizedHessianF,Inf))
    end
    
    for k = 1:S.maxIter
        S.log("========================================================")
        S.log("Main path following scheme - Step $k")
        S.log("========================================================")
        
        try
            t += S.gamma / starnorm(S.c, I, S.solveLS)
        catch e
            if isa(e, SingularException)
                S.log("Hessian singular at t = $t < $(S.tolInv)")
                I.msg *= " Hessian singular  - "
                iterationCount = -k
                S.eps = S.tolFactor / t    
                break
            else
                rethrow()
            end
        end
        
        I.dx = S.solveLS(I.hessianF, t * S.c + I.gradientF, I.P)
        apply_descent!(I, S, useBacktracking=false)
        
        if S.logObjective
            write_step(S.mainObjectiveFileName, k, S.tolFactor / t)
        end

        if S.logCondition
            write_step(S.mainConditionFileName, k, cond(I.unfactorizedHessianF,Inf))
        end
        
        if t >= S.tolInv
            S.log("t = $t ≥ $(S.tolInv)")
            iterationCount = k
            break
        end
        S.log("t = $t < $(S.tolInv)")

        if k == S.maxIter
            S.log("Maximum number of iterations reached.")
            I.msg *= " Exceeded iterations  - "
            iterationCount = -k
            break
        end
    end

    return iterationCount
end

"""
    select_pathfollowing(stepsize::Stepsize) -> Function, Function
    
Returns functions for auxilliary and main pathfollowing dependend on the stepsize.
"""
function select_pathfollowing(stepsize::Stepsize)  
    if stepsize === LONG
        return pathfollowing_auxilliary_long, pathfollowing_main_long
    elseif stepsize === ADAPTIVE
        return pathfollowing_auxilliary_adaptive, pathfollowing_main_adaptive
    else
        return pathfollowing_auxilliary_short, pathfollowing_main_short
    end
end
