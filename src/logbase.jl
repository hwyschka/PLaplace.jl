"""
$(TYPEDSIGNATURES)

Reimplemtation of deprecated function cpad based on current implementation of lpad and rpad.
"""
function cpad(
    s::Union{AbstractChar,AbstractString},
    n::Integer,
    p::Union{AbstractChar,AbstractString}=' ',
) :: String
    n = Int(n)::Int
    m = signed(n) - Int(textwidth(s))::Int
    m ≤ 0 && return string(s)
    
    l = textwidth(p)
    q, r = divrem(m, l)

    t, v = divrem(q, 2)
    lp = p^(trunc(Int, t))
    rp = p^(trunc(Int, t)+v)

    r == 0 ? string(lp, s, rp) : string(lp, s, rp, first(p, r))
end

cpad(s, n::Integer, p::Union{AbstractChar,AbstractString}=' ') = 
    cpad(string(s)::AbstractString, n, string(p))

"""
$(TYPEDSIGNATURES)

Prints a vertical line of width corresponding to the log to the given output.
"""
function log_print_line(out::Function)
    out("⧆" * repeat("=", 106) * "⧆\n")
end

log_print_footer(out::Function) = log_print_line(out)

"""
$(TYPEDSIGNATURES)

Prints a text shorter than the width of the log
right aligned with frame to the given output.
"""
function log_print_text(out::Function, text::String)
    out("‖" * rpad(text,106)* "‖\n")
end

"""
$(TYPEDSIGNATURES)

Prints a text shorter than the width of the log
center aligned with frame to the given output.
"""
function log_print_text_center(out::Function, text::String)
    out("‖" * cpad(text,106)* "‖\n")
end

"""
$(TYPEDSIGNATURES)

Prints spacing between log blocks to the given output.
"""
function log_print_spacing(out::Function)
    out("\n\n")
end

"""
$(TYPEDSIGNATURES)

Prints start of log to the given output.
"""
function log_print_start(out::Function)
    log_print_line(out)
    out("‖" * cpad("PLaplace.jl",106) * "‖\n")
    log_print_line(out)
end

"""
$(TYPEDSIGNATURES)
    
Prints problem data that had to be specified for the algorithm to the given output. 
"""
function log_print_defaultdata(
    out::Function,
    mesh::Mesh,
    eps::Float64,
    p::Float64,
    stepsize::String
)
    seps = @sprintf("%.5e", eps)
    sp = @sprintf("%.5f", p)

    c1 = 20
    c2 = 106 - 1 - c1
    
    log_print_spacing(out)
    log_print_line(out)

    out("‖" * rpad("Mesh - #nodes",c1))
    out("|" * rpad(string(mesh.nnodes),c2) * "‖\n")

    out("‖" * rpad("     - #elements",c1))
    out("|" * rpad(string(mesh.nelems),c2) * "‖\n")

    out("‖" * rpad("Accuracy ε",c1))
    out("|" * rpad(seps,c2) * "‖\n")

    out("‖" * rpad("PDE Parameter p",c1))
    out("|" * rpad(sp,c2) * "‖\n")

    out("‖" * rpad("Step Size",c1))
    out("|" * rpad(stepsize,c2) * "‖\n")

    log_print_line(out)
end

"""
$(TYPEDSIGNATURES)
    
Prints block header for the log to the given output. 
"""
function log_print_header(out::Function)
    log_print_spacing(out)
    log_print_line(out)

    out("‖" * cpad("Iteration",17))
    out("‖" * cpad("Step",27))
    out("‖" * cpad("Backtracking",12))
    out("‖" * cpad("Criteria",47) * "‖\n")


    out("‖" * cpad("Total",8))
    out("|" * cpad("Section",8))

    out("‖" * cpad("??",2))
    out("|" * cpad("LA",2))
    out("|" * cpad("κ",9))
    out("|" * cpad("‖⋅‖ !< β",11))

    out("‖" * cpad("it",2))
    out("|" * cpad("r",9))

    out("‖" * cpad("t",15))
    out("|" * cpad("‖F'(x)‖",15))
    out("|" * cpad("Crit Bound",15) * "‖\n")

    log_print_line(out)
end

"""
$(TYPEDSIGNATURES)
    
Prints one iteration for the log to the given output.
In particular contains all relevant information and format parsing.
"""
function log_print_iteration(
    out::Function,
    itr::Union{Int64,Missing},
    offset::Int64,
    prefix::String,
    type::Union{String,Missing},
    lastaccept::Union{Int64,Missing},
    kappa::Union{Float64,Missing},
    accnorm::Union{Float64,Missing},
    searchitr::Union{Int64,Missing},
    searchval::Union{Float64,Missing},
    crit::Union{Float64,Missing},
    t::Union{Float64,Missing},
    bound::Float64
)
    if !ismissing(itr) && !iszero(itr) && (itr % 10) == 0 
        log_print_footer(out)
        log_print_header(out)
    end

    sit = ismissing(itr) ? repeat("/",8) : lpad(string(itr+offset),8)
    si = ismissing(itr) ? repeat("/",8) : prefix * lpad(string(itr),7)

    sty = ismissing(type) ? repeat("/",2) : rpad(type,2)
    sla = ismissing(lastaccept) ? repeat("/",2) : lpad(string(lastaccept),2)
    sk = ismissing(kappa) ? repeat("/",9) : lpad(@sprintf("%.3e", kappa),9)
    san = ismissing(accnorm) ? repeat("/",11) : lpad(@sprintf("%.5e", accnorm),11)

    stri = ismissing(searchitr) ? repeat("/",2) : lpad(string(searchitr),2)
    strv = ismissing(searchval) ? repeat("/",9) : lpad(@sprintf("%.3e", searchval),9)

    st = ismissing(t) ? repeat("/",15) : @sprintf("%.9e", t)
    sc = ismissing(crit) ? repeat("/",15) : @sprintf("%.9e", crit)
    sb = @sprintf("%.9e", bound)


    out("‖" * sit)
    out("|" * si)
    
    out("‖" * sty)
    out("|" * sla)
    out("|" * sk)
    out("|" * san)
    
    out("‖" * stri)
    out("|" * strv)

    out("‖" * lpad(st,15))
    out("|" * lpad(sc,15))
    out("|" * lpad(sb,15) * "‖\n")
end

"""
$(TYPEDSIGNATURES)
    
Prints setup updated designed to fit the log to the given output.
"""
function log_print_setup(out::Function, stage::Int64)
    c1 = 40
    c2 = 106 - 1 - c1
    if stage == 0
        log_print_spacing(out)
        log_print_line(out)
        out("‖" * rpad("Evaluating Functions",c1))
    elseif stage == 1
        out("|" * rpad(" ✓",c2) * "‖\n")
    elseif stage == 2
        out("‖" * rpad("Initializing Barrier",c1))
    elseif stage == 3
        out("|" * rpad(" ✓",c2) * "‖\n")
    elseif stage == 4
        out("‖" * rpad("Assembling Static Data",c1))
    elseif stage == 5
        out("|" * rpad(" ✓",c2) * "‖\n")
    elseif stage == 6
        out("‖" * rpad("Initializing Iteration Data",c1))
    elseif stage == 7
        out("|" * rpad(" ✓",c2) * "‖\n")
        log_print_line(out)
    end
end

"""
$(TYPEDSIGNATURES)
    
Prints setup updated designed to fit the log to the given output.
"""
function log_print_change_factorization(out::Function)
    log_print_line(out)
    log_print_text(
        out,
        "Hessian not numerical sym. pos. def. Changed factorization to LU."
    )
    log_print_line(out)
end

"""
$(TYPEDSIGNATURES)
    
Prints setup updated designed to fit the log to the given output.
"""
function log_print_change_preconditioner(out::Function)
    log_print_line(out)
    log_print_text(
        out,
        "Hessian not numerical pos def. Changed precondtioner to LU."
    )
    log_print_line(out)
end

"""
$(TYPEDSIGNATURES)
    
Prints statistics of the full solution process to the given output.
In particular contains required iterations and time as well as termination messages.
"""
function log_print_statistics(
    out::Function,
    tsetup::Union{Float64,Missing},
    taux::Union{Float64,Missing},
    tmain::Union{Float64,Missing},
    Naux::Union{Int64,Missing},
    Nmain::Union{Int64,Missing},
    eps::Float64,
    msg::String
)
    stsetup = ismissing(tsetup) ? repeat("-",5) : @sprintf("%.3f", tsetup)
    staux = ismissing(taux) ? repeat("-",5) : @sprintf("%.3f", taux)
    stmain = ismissing(tmain) ? repeat("-",5) : @sprintf("%.3f", tmain)
    
    tsum = 0
    tsum += ismissing(tsetup) ? 0 : tsetup
    tsum += ismissing(taux) ? 0 : taux
    tsum += ismissing(tmain) ? 0 : tmain
    stsum = iszero(tsum) ? repeat("-",5) : @sprintf("%.3f", tsum)

    sNaux = ismissing(Naux) ? repeat("-",3) : string(abs(Naux))
    sNmain = ismissing(Nmain) ? repeat("-",3) : string(abs(Nmain))

    N = 0
    N += ismissing(Naux) ? 0 : abs(Naux)
    N += ismissing(Nmain) ? 0 : abs(Nmain)
    sN = iszero(N) ? repeat("-",4) : string(N)

    c1 = 15
    c2 = 12
    c3 = 106 - 2 - c1 - c2
    
    log_print_spacing(out)
    log_print_line(out)

    out("‖" * rpad("",c1))
    out("|" * rpad("Iterations",c2))
    out("|" * rpad("Time (s)",c3) * "‖\n")

    log_print_line(out)

    out("‖" * rpad("Setup",c1))
    out("|" * rpad("--",c2))
    out("|" * rpad(stsetup,c3) * "‖\n")

    out("‖" * rpad("Auxiliary",c1))
    out("|" * rpad(sNaux,c2))
    out("|" * rpad(staux,c3) * "‖\n")

    out("‖" * rpad("Main",c1))
    out("|" * rpad(sNmain,c2))
    out("|" * rpad(stmain,c3) * "‖\n")

    out("‖" * rpad("Combined",c1))
    out("|" * rpad(sN,c2))
    out("|" * rpad(stsum,c3) * "‖\n")

    out("‖" * rpad("Message",c1))
    out("|" * rpad(msg,c2+1+c3)* "‖\n")

    log_print_line(out)
end

"""
$(TYPEDSIGNATURES)
    
Clears and writes a header for an stepwise log to the given file. 
"""
function log_write_header(filename::String)
    open(filename, "w") do file
        write(file, rpad("Iteration",18), "|")
        write(file, rpad("Step",27), "|")
        write(file, rpad("Backtracking",12), "|")
        write(file, rpad("Criteria",48), "|")
        write(file, rpad("Hessian",15), "\n")


        write(file, rpad("Total",8), "|")
        write(file, rpad("P",1), "|")
        write(file, rpad("Part",7), "|")

        write(file, rpad("St",2), "|")
        write(file, rpad("LA",2), "|")
        write(file, rpad("kappa",9), "|")
        write(file, rpad("centering",11), "|")

        write(file, rpad("it",2), "|")
        write(file, rpad("r",9), "|")

        write(file, rpad("t",16), "|")
        write(file, rpad("Grad Norm",15), "|")
        write(file, rpad("Crit Bound",15), "|")

        write(file, rpad("Condition",15), "\n")

        write(file, repeat("-", 124), "\n")
        write(file, "\$Iterations", "\n")
    end
end

"""
$(TYPEDSIGNATURES)
    
Writes one iteration for an stepwise log to the given file.
In particular contains all relevant information and format parsing.
"""
function log_write_iteration(
    filename::String,
    itr::Union{Int64,Missing},
    offset::Int64,
    prefix::String,
    type::Union{String,Missing},
    lastaccept::Union{Int64,Missing},
    kappa::Union{Float64,Missing},
    accnorm::Union{Float64,Missing},
    searchitr::Union{Int64,Missing},
    searchval::Union{Float64,Missing},
    crit::Union{Float64,Missing},
    t::Union{Float64,Missing},
    bound::Float64,
    condition::Union{Float64,Missing}
)
    sitrtot = ismissing(itr) ? repeat("~",8) : @sprintf("%.08i", itr+offset)
    spre = prefix 
    sitr = ismissing(itr) ? repeat("~",7) : @sprintf("%.07i", itr)

    sty = ismissing(type) ? repeat("~",2) : rpad(type, 2, '~')
    sla = ismissing(lastaccept) ? repeat("~",2) : @sprintf("%.02i", lastaccept)
    sk = ismissing(kappa) ? repeat("~",9) : @sprintf("%.3e", kappa)
    san = ismissing(accnorm) ? repeat("~",11) : @sprintf("%.5e", accnorm)

    stri = ismissing(searchitr) ? repeat("~",2) : @sprintf("%.02i", searchitr)
    strv = ismissing(searchval) ? repeat("~",9) : @sprintf("%.3e", searchval)

    st = ismissing(t) ? repeat("~",16) : @sprintf("%+.9e", t)
    sc = ismissing(crit) ? repeat("~",15) : @sprintf("%.9e", crit)
    sb = @sprintf("%.9e", bound)

    scond = ismissing(condition) ? repeat("~",15) : @sprintf("%.9e", condition)


    open(filename, "a") do file
        write(file, sitrtot, " ")
        write(file, spre, " ")
        write(file, sitr, " ")

        write(file, sty, " ")
        write(file, sla, " ")
        write(file, sk, " ")
        write(file, san, " ")

        write(file, stri, " ")
        write(file, strv, " ")

        write(file, st, " ")
        write(file, sc, " ")
        write(file, sb, " ")

        write(file, scond, "\n")
    end
end
