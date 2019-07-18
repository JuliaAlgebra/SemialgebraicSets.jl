export @set

struct PolynomialInequality{PT<:AbstractPolynomialLike}
    p::PT
end
function Base.intersect(ineq::PolynomialInequality; kws...)
    return BasicSemialgebraicSet(FullSpace(), [ineq.p])
end

struct PolynomialEquality{PT<:AbstractPolynomialLike}
    p::PT
end
function equality(lhs, rhs)
    return PolynomialEquality(lhs - rhs)
end
function Base.intersect(eq::PolynomialEquality; lib_or_solver=nothing)
    if lib_or_solver === nothing
        return algebraicset([eq.p])
    else
        return algebraicset([eq.p], lib_or_solver)
    end
end

const Element = Union{PolynomialInequality, PolynomialEquality, FixedVariable}

function Base.intersect(el::Element,
                        args...; kws...)
    return intersect(intersect(el; kws...), args...; kws...)
end

function Base.intersect(set::AbstractBasicSemialgebraicSet,
                        el::Element, args...; kws...)
    return intersect(set, intersect(el; kws...), args...; kws...)
end

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

function appendconstraints!(elements, expr, _error)
    if Base.Meta.isexpr(expr, :call)
        try
            sense, vectorized = _canonicalize_sense(expr.args[1], _error)
            @assert !vectorized
            if sense == :(>=)
                push!(elements, esc(:(SemialgebraicSets.PolynomialInequality($(expr.args[2]) - $(expr.args[3])))))
            elseif sense == :(<=)
                push!(elements, esc(:(SemialgebraicSets.PolynomialInequality($(expr.args[3]) - $(expr.args[2])))))
            elseif sense == :(==)
                push!(elements, esc(:(SemialgebraicSets.equality($(expr.args[2]), $(expr.args[3])))))
            else
                _error("Unrecognized sense $(string(sense)) in domain specification")
            end
        catch
            push!(elements, esc(expr))
        end
    elseif Base.Meta.isexpr(expr, :&&)
        map(t -> appendconstraints!(elements, t, _error), expr.args)
    else
        push!(elements, esc(expr))
    end
    return nothing
end

macro set(expr, library=nothing)
    elements = []
    appendconstraints!(elements, expr, msg -> error("In @set($expr: ", msg))
    if library === nothing
        return :( intersect($(elements...)) )
    else
        return :( intersect($(elements...); lib_or_solver=$(esc(library))) )
    end
end
