"""
$(TYPEDEF)

Structure to hold general information considering the log.
In particular if the algorithm is verbose and if a log file is generated.
    
# Fields
$(TYPEDFIELDS)
"""
mutable struct LogData
    "Flag indicating if output is verbose."
    isverbose::Bool
    
    "Output function is print() if verbose otherwise does nothing."
    out::Function
    
    "Flag for logging the objective to file during the iteration."
    logobjective::Bool
    
    "Flag for logging the condition of the linear system to file during the iteration."
    trackcondition::Bool
    
    "File name for writing log to file."
    file::Union{String,Missing}
end

"""
$(TYPEDSIGNATURES)

Constructor for [LogData](@ref) with default values.
"""
function LogData()
    return LogData(true, print, false, false, missing)
end

"""
$(TYPEDSIGNATURES)

Constructor for [LogData](@ref) with data.
"""
function LogData(
    verbose::Bool,
    filename::Union{String,Missing},
    logcondition::Bool
) 
    logobjective = !ismissing(filename)
    
    fn::Union{String,Missing} = missing
    if logobjective
        fn = endswith(".txt", filename) ? filename : filename * ".txt"

        log_write_header(fn)
    end

    trackcondition = logobjective && logcondition

    return LogData(
        verbose,
        verbose ? print : emptyfunction,
        logobjective,
        trackcondition,
        fn
    )
end

"""
$(TYPEDSIGNATURES)

Handles horizontal line in the log for output stream. Is not relevant for file.
"""
function log_line(data::LogData)
    data.isverbose && log_print_line(data.out)
end

"""
$(TYPEDSIGNATURES)

Handles footer in the log for output stream. Is not relevant for file.
"""
function log_footer(data::LogData)
    data.isverbose && log_print_footer(data.out)
end

"""
$(TYPEDSIGNATURES)

Handles text line in the log for output stream. Is not relevant for file.
"""
function log_text(data::LogData, msg::String)
    data.isverbose && log_print_text(data.out, msg)
end

"""
$(TYPEDSIGNATURES)

Handles start block in the log for output stream. Is not relevant for file.
"""
function log_start(data::LogData)
    data.isverbose && log_print_start(data.out)
end

"""
$(TYPEDSIGNATURES)
    
Handles default data in the log for output stream. Is not relevant for file.
"""
function log_defaultdata(
    data::LogData,
    mesh::Mesh,
    eps::Float64,
    p::Float64,
    stepsize::String
)
    if data.isverbose
        log_print_defaultdata(data.out, mesh, eps, p, stepsize)
    end
end

"""
$(TYPEDSIGNATURES)
    
Handles header in the log for output stream. Is not relevant for file.
"""
function log_header(data::LogData)
    data.isverbose && log_print_header(data.out)
end

log_inital(data::LogData) = log_header(data)


"""
$(TYPEDSIGNATURES)
    
Handles an iteration in the log for output stream and file generation.
"""
function log_iteration(
    data::LogData,
    itr::Union{Int64,Missing},
    offset::Int64,
    prefix::String,
    type::Union{String,Missing},
    lastaccept::Union{Int64,Missing},
    kappa::Union{Float64,Missing},
    accnorm::Union{Float64,Missing},
    searchitr::Union{Int64,Missing},
    searchval::Union{Float64,Missing},
    condition::Union{Float64,Missing},
    crit::Union{Float64,Missing},
    t::Union{Float64,Missing},
    bound::Float64
)
    if data.logobjective
        log_write_iteration(
            data.file,
            itr, offset, prefix,
            type, lastaccept, kappa, accnorm,
            searchitr, searchval,
            crit, t, bound,
            condition
        )
    end

    if data.isverbose
        log_print_iteration(
            data.out,
            itr,offset,prefix,
            type,lastaccept,kappa,accnorm,
            searchitr,searchval,
            crit,t,bound
        )
    end
end

"""
$(TYPEDSIGNATURES)
    
Handles setup in the log for output stream. Is not relevant for file.
"""
function log_setup(data::LogData, stage::Int64)
    if data.isverbose 
        log_print_setup(data.out, stage)
    end
end

"""
$(TYPEDSIGNATURES)
    
Handles change of factorization in the log for output stream. Is not relevant for file.
"""
function log_change_factorization(data::LogData)
    if data.isverbose 
        log_print_change_factorization(data.out)
    end
end

"""
$(TYPEDSIGNATURES)
    
Handles change of preconditioner in the log for output stream. Is not relevant for file.
"""
function log_change_preconditioner(data::LogData)
    if data.isverbose
        log_print_change_preconditioner(data.out)
    end
end

"""
$(TYPEDSIGNATURES)
    
Handles statistics in the log for output stream. Is not relevant for file.
"""
function log_statistics(
    data::LogData,
    tsetup::Union{Float64,Missing},
    taux::Union{Float64,Missing},
    tmain::Union{Float64,Missing},
    Naux::Union{Int64,Missing},
    Nmain::Union{Int64,Missing},
    eps::Float64,
    msg::String
)
    if data.isverbose
        log_print_statistics(
            data.out,
            tsetup, taux, tmain,
            Naux, Nmain,
            eps,
            msg
        )
    end
end
