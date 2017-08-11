export @set

# Taken from JuMP/macros.jl
function _canonicalize_sense(sns::Symbol, _error)
    if sns == :(==)
        return (:(==),false)
    elseif sns == :(>=) || sns == :(≥)
        return (:(>=),false)
    elseif sns == :(<=) || sns == :(≤)
        return (:(<=),false)
    elseif sns == :(.==)
        return (:(==),true)
    elseif sns == :(.>=) || sns == :(.≥)
        return (:(>=),true)
    elseif sns == :(.<=) || sns == :(.≤)
        return (:(<=),true)
    else
        _error("Unrecognized sense $sns")
    end
end

function appendconstraints!(domains, domaineqs, domainineqs, expr, _error)
    if Base.Meta.isexpr(expr, :call)
        try
            sense, vectorized = _canonicalize_sense(expr.args[1], _error)
            @assert !vectorized
            if sense == :(>=)
                push!(domainineqs, esc(:($(expr.args[2]) - $(expr.args[3]))))
            elseif sense == :(<=)
                push!(domainineqs, esc(:($(expr.args[3]) - $(expr.args[2]))))
            elseif sense == :(==)
                push!(domaineqs, esc(:($(expr.args[2]) - $(expr.args[3]))))
            else
                _error("Unrecognized sense $(string(sense)) in domain specification")
            end
        catch
            push!(domains, esc(expr))
        end
    elseif Base.Meta.isexpr(expr, :&&)
        map(t -> appendconstraints!(domains, domaineqs, domainineqs, t, _error), expr.args)
    else
        push!(domains, esc(expr))
    end
    nothing
end

function builddomain(domains, domaineqs, domainineqs)
    domainaffs = gensym()

    PT = gensym()
    T = gensym()
    if isempty(domaineqs) && isempty(domainineqs)
        if isempty(domains)
            code = :( $domainaffs = FullSpace() )
        elseif length(domains) == 1
            code = :( $domainaffs = $(domains[1]) )
        else
            code = :( $domainaffs = intersect($(domains...)) )
        end
    else
        eqs = gensym()
        ineqs = gensym()
        code = quote
            $eqs = tuple($(domaineqs...))
            $ineqs = tuple($(domainineqs...))
            $PT = Base.promote_typeof($eqs..., $ineqs...)
            $T = coefficienttype($PT)
        end

        lin = gensym()
        code = :( $code; $lin = AlgebraicSet{$T, $PT}($PT[$eqs...]) )
        basic = gensym()
        if isempty(domainineqs)
            code = :( $code; $basic = $lin )
        else
            code = :( $code; $basic = BasicSemialgebraicSet{$T, $PT}($lin, $PT[$ineqs...]) )
        end

        if !isempty(domains)
            code = :( $code; $domainaffs = intersect($basic, $(domains...)) )
        else
            code = :( $code; $domainaffs = $basic )
        end
    end
    domainaffs, code
end

macro set(expr)
    domains = []
    domaineqs = []
    domainineqs = []
    appendconstraints!(domains, domaineqs, domainineqs, expr, msg -> error("In @set($expr: ", msg))
    domainvar, domaincode = builddomain(domains, domaineqs, domainineqs)
    quote
        $domaincode
        $domainvar
    end
end
