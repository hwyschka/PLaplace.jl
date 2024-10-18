function isadmissible_finite(x::AbstractVector{Float64}, S::StaticData)
    z = zeros(Float64, S.m)
    for (key, val) in S.D
        z -= (val * x[1:S.lengthu] + S.b[key]).^2
    end

    s = x[(S.lengthu + 1):(S.lengthu + S.m)]

    for si in s
        si <= 0 && return false
    end
        
    tau = S.R .- (S.omega .* s)
    for taui in tau
        taui <= 0 && return false
    end

    z += s.^(2 / S.p)
    for zi in z
        zi <= 0 && return false
    end
    
    return true
end

function compute_value_finite(x::AbstractVector{Float64}, S::StaticData)
    s = x[(S.lengthu + 1):(S.lengthu + S.m)]
    tau = S.R .- (S.omega .* s)

    y = Dict{Tuple{Int64,Int64},AbstractVector{Float64}}()
    for (key, val) in S.D
        y[key] = val * x[1:S.lengthu] + S.b[key]
    end

    z = zeros(Float64, S.m)
    z += s.^(2 / S.p)
    for (key, val) in y
        z -= val.^2
    end

    return -sum(log.(z)) - sum(log.(tau)) - S.alpha*sum(log.(s))
end

function compute_gradient_finite(x::AbstractVector{Float64}, S::StaticData)
    s = x[(S.lengthu + 1):(S.lengthu + S.m)]
    tau = S.R .- (S.omega .* s)

    y = Dict{Tuple{Int64,Int64},AbstractVector{Float64}}()
    for (key, val) in S.D
        y[key] = val * x[1:S.lengthu] + S.b[key]
    end

    z = zeros(Float64, S.m)
    z += s.^(2 / S.p)
    for (key, val) in y
        z -= val.^2
    end

    Fu = zeros(Float64, S.lengthu)
    for (key, val) in S.D
        Fu += val' * (y[key] ./ z)
    end
    Fu *= 2

    Fs = zeros(Float64, S.m)
    Fs -= 2 / S.p * z.^(-1) .* s.^(2.0 / S.p - 1) 
    Fs += S.omega ./ tau 
    Fs -= S.alpha ./ s

    return [Fu; Fs]
end

function compute_hessian_finite(x::AbstractVector{Float64}, S::StaticData)
    s = x[(S.lengthu + 1):(S.lengthu + S.m)]
    tau = S.R .- (S.omega .* s)

    y = Dict{Tuple{Int64,Int64},AbstractVector{Float64}}()
    for (key, val) in S.D
        y[key] = val * x[1:S.lengthu] + S.b[key]
    end    

    z = zeros(Float64, S.m)
    z += s.^(2 / S.p)
    for (key, val) in y
        z -= val.^2
    end

    F_uu1 = spzeros(Float64,S.lengthu,S.lengthu)
    F_uu2 = spzeros(Float64,S.lengthu,S.lengthu)
    for (key1, val1) in S.D
        F_uu1 += val1' * Diagonal(z.^(-1)) * val1
        for (key2, val2) in S.D
            F_uu2 += (Diagonal(y[key1]) * val1)' * Diagonal(z.^(-2)) * (Diagonal(y[key2]) * val2)
        end
    end
    F_uu = 2 * F_uu1 + 4 * F_uu2
    
    Fus = spzeros(S.lengthu, S.m)
    for (key, val) in S.D
        Fus -= (Diagonal(y[key]) * val)' * Diagonal(z.^(-2)) * Diagonal(s.^(2.0/S.p - 1))
    end
    Fus *= 4 / S.p

    F_ss = zeros(S.m)
    F_ss -= 2 / S.p * (2.0 / S.p - 1.0) .* z.^(-1) .* s.^(2.0 / S.p - 2)
    F_ss += 4 / S.p^2 .* z.^(-2) .* s.^(4.0 / S.p - 2)
    F_ss += S.alpha * s.^(-2)
    F_ss += S.omega.^2 .* tau.^(-2)
    
    hessF = [Symmetric(F_uu) Fus; Fus' Diagonal(F_ss)]

    return hessF
end

function compute_initialguess_finite(S::StaticData)
    t = zeros(Float64, S.m)
    for (key, val) in S.b
        t += val.^2
    end

    s = 1 .+ t.^(S.p/2)

    return [zeros(Float64, S.lengthu); s]
end
