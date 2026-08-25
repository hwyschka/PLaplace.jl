"""
$(TYPEDEF)

Structure to store statistics about a PLaplace run.
Supposed to be generated as a lighter version of PLaplaceData,
in particular does not carry problem, mesh and solution.
Is also used to post-process statistic logs.

# Fields
$(TYPEDFIELDS)
"""
mutable struct PLaplaceRunStatistics
    "PDE parameter."
    p::Float64

    "Accuracy of the solution in a variational sense."
    eps::Float64

    "Accuracy of the solution in a variational sense."
    n::Int64

    "Accuracy of the solution in a variational sense."
    m::Int64

    "Stepping scheme of the interior point method."
    stepsize::Stepsize
    
    "Number of iterations in the auxiliary path-following.
        Is missing if run did not reach auxiliary stage."
    Naux::Union{Int64,Missing}

    "Number of iterations for the main path-following.
        Is missing if run did not reach main stage."
    Nmain::Union{Int64,Missing}

    "Number of iterations for the auxiliary and main path-following combined.
        Is missing if run did not reach auxiliary stage."
    Nsum::Union{Int64,Missing}

    "Time required for the setup.
        Is missing if run did not reach setup stage."
    tsetup::Union{Float64,Missing}

    "Time required for the auxiliary path-following.
        Is missing if run did not reach auxiliary stage."
    taux::Union{Float64,Missing}

    "Time required for the main path-following.
        Is missing if run did not reach main stage."
    tmain::Union{Float64,Missing}

    "Time required for the setup as well as auxiliary and main path-following combined.
        Is missing if run did not reach setup."
    tsum::Union{Float64,Missing}

    "Notifications from the iteration. In particular contains changes of solvers and
        preconditioners and information on early stops."
    msg::String
end

function PLaplaceRunStatistics(data::PLaplaceData)
    Nsum::Union{Int64, Missing} = missing
    if !ismissing(data.Naux)
        Nsum = data.Naux
        if !ismissing(data.Nmain)
            Nsum += abs(data.Nmain)
            Nsum *= sign(data.Nmain)
        end
    end

    tsum::Union{Float64, Missing} = missing
    if !ismissing(data.tsetup)
        tsum = data.tsetup
        tsum += ismissing(data.taux) ? 0 : data.taux
        tsum += ismissing(data.tmain) ? 0 : data.tmain
    end

    return PLaplaceRunStatistics(
        data.p,
        data.eps,
        data.mesh.nnodes,
        data.mesh.nelems,
        data.stepsize,
        data.Naux,
        data.Nmain,
        Nsum,
        data.tsetup,
        data.taux,
        data.tmain,
        tsum,
        data.msg
    )
end

"""
    write_statistics_header(filename::String; guarded::Bool=false)
    
Clears and writes a header for a statistics log to the given file.
In case the file does not exist, it will be created, but only if the path exists.
If guarded checks before if file already contains a header
and then does not overwrite potentially previous results.
"""
function write_statistics_header(filename::String; guarded::Bool=false)
    fn = occursin(".", filename) ? filename : filename * ".txt"

    if guarded
        check_statistics_header(fn) && return
    end

    open(fn, "w") do file
        write(file, rpad("p",6), "|")
        write(file, rpad("eps",13), "|")
        write(file, rpad("n",7), "|")
        write(file, rpad("m",7), "|")
        write(file, rpad("Scheme",8), "|")
        write(file, rpad("Naux",6), "|")
        write(file, rpad("Nmain",6), "|")
        write(file, rpad("N",7), "|")
        write(file, rpad("t setup",9), "|")
        write(file, rpad("t aux",9), "|")
        write(file, rpad("t main",9), "|")
        write(file, rpad("t sum",9), "|")
        write(file, rpad("Message",10), "\n")
        write(file, repeat("-", 115), "\n")
        write(file, "\$Simulations", "\n")
    end
end

"""
$(TYPEDSIGNATURES)
    
Checks if given file exists and already contains a statistics header. 
"""
function check_statistics_header(filename::String) :: Bool
    fn = occursin(".", filename) ? filename : filename * ".txt"

    if !isfile(fn)
        return false
    end

    f = open(fn)
    l = readline(f)
    close(f)
    a = split(l, "|")
    
    !contains(a[1],"p") && return false
    !contains(a[2],"eps") && return false
    !contains(a[13],"Message") && return false

    return true
end

"""
$(TYPEDSIGNATURES)
    
Writes statistics line corresponding to the header to a given log file. 
"""
function write_statistics(filename::String, data::PLaplaceRunStatistics)
    fn = occursin(".", filename) ? filename : filename * ".txt"

    sp = @sprintf("%06.3f", data.p)
    se = @sprintf("%.7e", data.eps)
    sn = @sprintf("%.07i", data.n)
    sm = @sprintf("%.07i", data.m)

    st = replace(@sprintf("%-8s", data.stepsize)," "=>"~")
    na = ismissing(data.Naux) ? repeat("~", 6) : @sprintf("%+.05i", data.Naux)
    nm = ismissing(data.Nmain) ? repeat("~", 6) : @sprintf("%+.05i", data.Nmain)
    nc = ismissing(data.Nsum) ? repeat("~", 7) : @sprintf("%+.06i", data.Nsum)

    ts = ismissing(data.tsetup) ? "~~~~~.~~~" : @sprintf("%09.3f", data.tsetup)
    ta = ismissing(data.taux) ? "~~~~~.~~~" : @sprintf("%09.3f", data.taux)
    tm = ismissing(data.tmain) ? "~~~~~.~~~" : @sprintf("%09.3f", data.tmain)
    tc = ismissing(data.tsum) ? "~~~~~.~~~" : @sprintf("%09.3f", data.tsum)

    open(fn, "a") do file
        write(file, sp, " ")
        write(file, se, " ")
        write(file, sn, " ")
        write(file, sm, " ")
        write(file, st, " ")
        write(file, na, " ")
        write(file, nm, " ")
        write(file, nc, " ")
        write(file, ts, " ")
        write(file, ta, " ")
        write(file, tm, " ")
        write(file, tc, " ")
        write(file, data.msg, "\n")
    end
end

write_statistics(filename::String, data::PLaplaceData) = 
    write_statistics(filename, PLaplaceRunStatistics(data))

"""
$(TYPEDSIGNATURES)
    
Reads given statistics file and returns an array of all logged runs. 
"""
function read_statistics(filename::String)
    
    list = Array{PLaplaceRunStatistics,1}()

    f = open(filename)

    while (!eof(f) && (l = readline(f)) != "\$Simulations")
    end

    while (!eof(f))
        l = readline(f)
        a = split(l, " ")

        msg = ""
        for i = 13:length(a)
            msg *= "$(a[i]) "    
        end

        run = PLaplaceRunStatistics(
            occursin("Inf", a[1]) ? Inf : parse(Float64, a[1]),
            parse(Float64, a[2]),
            parse(Int64, a[3]),
            parse(Int64, a[4]),
            parse(Stepsize, replace(a[5],"~"=>"")),
            
            occursin("~", a[6]) ? missing : parse(Int64, a[6]),
            occursin("~", a[7]) ? missing : parse(Int64, a[7]),
            occursin("~", a[8]) ? missing : parse(Int64, a[8]),

            occursin("~", a[9]) ? missing : parse(Float64, a[9]),
            occursin("~", a[10]) ? missing : parse(Float64, a[10]),
            occursin("~", a[11]) ? missing : parse(Float64, a[11]),
            occursin("~", a[12]) ? missing : parse(Float64, a[12]),
            msg
        )

        push!(list, run)
    end

    close(f)

    return list
end
