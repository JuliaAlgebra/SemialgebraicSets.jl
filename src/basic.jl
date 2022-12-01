export inequalities, basic_semialgebraic_set

struct BasicSemialgebraicSet{T,PT<:APL{T},AT<:AbstractAlgebraicSet} <:
       AbstractBasicSemialgebraicSet
    V::AT
    p::Vector{PT}
end
function BasicSemialgebraicSet{T,PT}() where {T,PT<:APL{T}}
    return BasicSemialgebraicSet(AlgebraicSet{T,PT}(), PT[])
end
function BasicSemialgebraicSet(
    V::AlgebraicSet{T,PT,A,S},
    p::Vector{PT},
) where {T,PT<:APL{T},A,S<:AbstractAlgebraicSolver}
    return BasicSemialgebraicSet{T,PT,typeof(V)}(V, p)
end
function BasicSemialgebraicSet(
    V::AlgebraicSet{T,PT,A,SO,U},
    p::Vector{PS},
) where {T,PT<:APL{T},S,PS<:APL{S},A,SO<:AbstractAlgebraicSolver,U}
    ST = promote_type(T, S)
    PST = promote_type(PT, PS)
    return BasicSemialgebraicSet(
        convert(AlgebraicSet{ST,PST,A,SO,U}, V),
        Vector{PST}(p),
    )
end
#BasicSemialgebraicSet{T, PT<:APL{T}}(V::AlgebraicSet{T, PT}, p::Vector{PT}) = BasicSemialgebraicSet{T, PT}(V, p)
function basic_semialgebraic_set(V, p)
    return BasicSemialgebraicSet(V, p)
end

function MP.changecoefficienttype(
    ::Type{BasicSemialgebraicSet{S,PS,AT}},
    T::Type,
) where {S,PS,AT}
    return BasicSemialgebraicSet{
        T,
        MP.changecoefficienttype(PS, T),
        MP.changecoefficienttype(AT, T),
    }
end

function Base.convert(
    ::Type{BasicSemialgebraicSet{T,PT,AT}},
    set::BasicSemialgebraicSet,
) where {T,PT,AT}
    return BasicSemialgebraicSet{T,PT,AT}(set.V, set.p)
end

function MP.variables(
    S::BasicSemialgebraicSet{T,PT,FullSpace},
) where {T,PT<:APL{T}}
    return MP.variables(S.p)
end
function MP.variables(S::BasicSemialgebraicSet)
    return sort(union(MP.variables(S.V), MP.variables(S.p)); rev = true)
end
nequalities(S::BasicSemialgebraicSet) = nequalities(S.V)
equalities(S::BasicSemialgebraicSet) = equalities(S.V)
add_equality!(S::BasicSemialgebraicSet, p) = add_equality!(S.V, p)
ninequalities(S::BasicSemialgebraicSet) = length(S.p)
inequalities(S::BasicSemialgebraicSet) = S.p
add_inequality!(S::BasicSemialgebraicSet, p) = push!(S.p, p)

function Base.intersect(S::BasicSemialgebraicSet, T::BasicSemialgebraicSet)
    return BasicSemialgebraicSet(S.V ∩ T.V, [S.p; T.p])
end
function Base.intersect(S::BasicSemialgebraicSet, T::AbstractAlgebraicSet)
    return BasicSemialgebraicSet(S.V ∩ T, copy(S.p))
end
Base.intersect(S::BasicSemialgebraicSet, ::FullSpace) = S
function Base.intersect(T::AbstractAlgebraicSet, S::BasicSemialgebraicSet)
    return intersect(S, T)
end

function Base.show(io::IO, V::BasicSemialgebraicSet)
    print(
        io,
        "{ (",
        join(variables(V), ", "),
        ") | ",
        join(string.(equalities(V)) .* " = 0", ", "),
    )
    if nequalities(V) > 0
        print(io, ", ")
    end
    return print(io, join(string.(inequalities(V)) .* " ≥ 0", ", "), " }")
end
function Base.show(io::IO, mime::MIME"text/plain", V::BasicSemialgebraicSet)
    print(io, "Basic semialgebraic Set defined by ")
    _show_els(io, "equalit", nequalities(V), equalities(V), "=")
    return _show_els(io, "inequalit", ninequalities(V), inequalities(V), "≥")
end
