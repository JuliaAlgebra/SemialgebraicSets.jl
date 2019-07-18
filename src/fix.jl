export FixedVariablesSet

struct FixedVariablesIdeal{V<:AbstractVariable, T<:Number, MT<:AbstractMonomialLike} <: AbstractPolynomialIdeal
    substitutions::Union{Nothing, Dict{V, T}}
end
# In that case, the ideal can generate any polynomial.
function generate_nonzero_constant(I::FixedVariablesIdeal)
    return I.substitutions === nothing
end
function Base.rem(p::AbstractPolynomialLike, I::FixedVariablesIdeal)
    if generate_nonzero_constant(I)
        return zero(p)
    else
        substitutions = [key => value for (key, value) in I.substitutions]
        return subs(p, substitutions...)
    end
end

struct FixedVariablesSet{V, T, MT} <: AbstractAlgebraicSet
    ideal::FixedVariablesIdeal{V, T, MT}
end
ideal(set::FixedVariablesSet, args...) = set.ideal
function nequalities(set::FixedVariablesSet)
    if set.ideal.substitutions === nothing
        return 1
    else
        return length(set.ideal.substitutions)
    end
end
function equalities(set::FixedVariablesSet{V, T, MT}) where {V, T, MT}
    if set.ideal.substitutions === nothing
        return [constantterm(one(T), MT)]
    else
        return [key - value for (key, value) in set.ideal.substitutions]
    end
end
function Base.intersect(V::FixedVariablesSet{V1, T1, MT1}, W::FixedVariablesSet{V2, T2, MT2}) where {V1, V2, T1, T2, MT1, MT2}
    # For `DynamicPolynomials`, they have the same type and for
    # `TypedPolynomials`, promoting would give `Monomial`.
    VT = V1 == V2 ? V1 : AbstractVariable
    T = promote_type(T1, T2)
    if V.ideal.substitutions === nothing || W.ideal.substitutions === nothing
        sub = nothing
    else
        has_dup = false
        function combine(a, b)
            if a != b
                has_dup = true
            end
            return a
        end
        sub = Dict{VT, T}()
        merge!(combine, sub, V.ideal.substitutions)
        merge!(combine, sub, W.ideal.substitutions)
        if has_dup
            sub = nothing
        end
    end
    ideal = FixedVariablesIdeal{VT, T, promote_type(MT1, MT2)}(sub)
    return FixedVariablesSet(ideal)
end
function Base.intersect(V::AlgebraicSet, W::FixedVariablesSet)
    return V ∩ algebraicset(W)
end
function Base.intersect(V::FixedVariablesSet, W::AlgebraicSet)
    # Need `W` to be first so that the intersection uses its solver.
    # and not the default one taken in `algebraicset(W)`.
    return W ∩ V
end

iszerodimensional(::FixedVariablesSet) = true
function Base.isempty(V::FixedVariablesSet)
    return generate_nonzero_constant(V.ideal)
end

# Assumes `isempty(V)`
function only_point(V::FixedVariablesSet)
    subs = collect(V.ideal.substitutions)
    sort!(subs, rev = true, by = sub -> sub[1])
    return [sub[2] for sub in subs]
end

function Base.iterate(V::FixedVariablesSet, state=nothing)
    if state === nothing && !isempty(V)
        return only_point(V), true
    else
        return nothing
    end
end
Base.length(V::FixedVariablesSet) = isempty(V) ? 0 : 1

struct FixedVariable{V<:AbstractVariable, T}
    variable::V
    value::T
end
function equality(variable::AbstractVariable, value::Number)
    return FixedVariable(variable, value)
end
function equality(value::Number, variable::AbstractVariable)
    return FixedVariable(variable, value)
end

function Base.intersect(el::FixedVariable; kws...)
    subs = Dict(el.variable => el.value)
    return FixedVariablesSet(FixedVariablesIdeal{
        typeof(el.variable), typeof(el.value), typeof(el.variable)}(subs))
end
