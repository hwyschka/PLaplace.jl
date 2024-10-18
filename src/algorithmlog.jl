"""
$(TYPEDEF)

Stores data read from an algorithm log file.
Can be used to generate plots from one solution process. 

# Fields
$(TYPEDFIELDS)
"""
mutable struct AlgorithmLogData
    "Running number of total iterations for each iteration."
    totalIterations::Array{Union{Int64,Missing},1}
    
    "Section tag in each iteration."
    sections::Array{String,1}
    
    "Running number of iterations in current section for each iteration."
    partIterations::Array{Union{Int64,Missing},1}

    "Tag for type of the update step for each iteration.
        Only relevant for long and adaptive stepping schemes."
    types::Array{Union{String,Missing},1}
    
    "Number of iterations since previous accept for each iteration
        Only relevant for adaptive stepping scheme."
    lastaccepts::Array{Union{Int64,Missing},1}
    
    "Parameter kappa for each iteration. Only relevant for adaptive stepping scheme."
    kappas::Array{Union{Float64,Missing},1}
    
    "Norm value for acceptance of step for each iteration.
        Only relevant for long and adaptive stepping schemes."
    accnorms::Array{Union{Float64,Missing},1}
    
    "Number of iterations in backtracking linesearch for each iteration.
        Only relevant for long and adaptive stepping schemes"
    searchIterations::Array{Union{Int64,Missing},1}
    
    "Scaling of the update generated by the backtracking for each iteration.
        Only relevant for long and adaptive stepping schemes."
    scalings::Array{Union{Float64,Missing},1}
    
    "Value of t for each iteration."
    ts::Array{Union{Float64,Missing},1}
    
    "Value of convergence criterum for each iteration."
    criteria::Array{Union{Float64,Missing},1}
    
    "Bound for convergence criteria for each iteration. Is constant per section."
    bounds::Array{Union{Float64,Missing},1}
    
    "Condition number of the system matrix for each iteration."
    conditions::Array{Union{Float64,Missing},1}
end

"""
    read_algorithmlog(filename::String)

Returns [AlgorithmLogData](@ref) by reading the given algorithm log file.
"""
function read_algorithmlog(filename::String)
    entries = 0
    open(filename) do f
        while (!eof(f) && (l = readline(f)) != "\$Iterations")
        end

        while (!eof(f))
            l = readline(f)
            entries += 1
        end
    end

    totalIterations = Array{Union{Int64,Missing},1}(undef,entries)
    sections = Array{String,1}(undef,entries)
    partIterations = Array{Union{Int64,Missing},1}(undef,entries)

    types = Array{Union{String,Missing},1}(undef,entries)
    lastaccepts = Array{Union{Int64,Missing},1}(undef,entries)
    kappas = Array{Union{Float64,Missing},1}(undef,entries)
    accnorms = Array{Union{Float64,Missing},1}(undef,entries)

    searchIterations = Array{Union{Int64,Missing},1}(undef,entries)
    scalings = Array{Union{Float64,Missing},1}(undef,entries)

    ts = Array{Union{Float64,Missing},1}(undef,entries)
    criteria = Array{Union{Float64,Missing},1}(undef,entries)
    bounds = Array{Union{Float64,Missing},1}(undef,entries)

    conditions = Array{Union{Float64,Missing},1}(undef,entries)

    f = open(fileName)

    while (!eof(f) && (l = readline(f)) != "\$Iterations")
    end

    k = 1
    while (!eof(f))
        l = readline(f)
        a = split(l, " ")

        totalIterations[k] = startswith(a[1], '~') ? missing : parse(Int64, a[1])
        sections[k] = a[2]
        partIterations[k] = startswith(a[3], '~') ? missing : parse(Int64, a[3])

        types[k] = startswith(a[4], '~') ? missing : a[4]
        lastaccepts[k] = startswith(a[5], '~') ? missing : parse(Int64, a[5])
        kappas[k] = startswith(a[6], '~') ? missing : parse(Float64, a[6])
        accnorms[k] = startswith(a[7], '~') ? missing : parse(Float64, a[7])

        searchIterations[k] = startswith(a[8], '~') ? missing : parse(Int64, a[8])
        scalings[k] = startswith(a[9], '~') ? missing : parse(Float64, a[9])

        ts[k] = startswith(a[10], '~') ? missing : parse(Float64, a[10])
        criteria[k] = startswith(a[11], '~') ? missing : parse(Float64, a[11])
        bounds[k] = startswith(a[12], '~') ? missing : parse(Float64, a[12])

        conditions[k] = startswith(a[13], '~') ? missing : parse(Float64, a[13])

        k += 1
    end

    return AlgorithmLogData(
        totalIterations,
        sections,
        partIterations,
        types,
        lastaccepts,
        kappas,
        accnorms,
        searchIterations,
        scalings,
        ts,
        criteria,
        bounds,
        conditions
    )
end
