export inequalities, basicsemialgebraicset

struct BasicSemialgebraicSet{T, PT<:APL{T}, AT <: AbstractAlgebraicSet} <: AbstractBasicSemialgebraicSet
    V::AT
    p::Vector{PT}
end
function BasicSemialgebraicSet{T, PT}() where {T, PT<:APL{T}}
    BasicSemialgebraicSet(AlgebraicSet{T, PT}(), PT[])
end
function BasicSemialgebraicSet(V::AlgebraicSet{T, PT, A, S}, p::Vector{PT}) where {T, PT<:APL{T}, A, S<:AbstractAlgebraicSolver}
    BasicSemialgebraicSet{T, PT, typeof(V)}(V, p)
end
function BasicSemialgebraicSet(V::AlgebraicSet{T, PT, A, ST}, p::Vector{PS}) where {T, PT<:APL{T}, S, PS<:APL{S}, A, ST<:AbstractAlgebraicSolver}
    U = promote_type(T, S)
    PU = promote_type(PT, PS)
    BasicSemialgebraicSet(convert(AlgebraicSet{U, PU, A, ST}, V), Vector{PU}(p))
end
#BasicSemialgebraicSet{T, PT<:APL{T}}(V::AlgebraicSet{T, PT}, p::Vector{PT}) = BasicSemialgebraicSet{T, PT}(V, p)
function basicsemialgebraicset(V, p)
    BasicSemialgebraicSet(V, p)
end

MP.changecoefficienttype(::Type{BasicSemialgebraicSet{S, PS, AT}}, T::Type) where {S, PS, AT} = BasicSemialgebraicSet{T, MP.changecoefficienttype(PS, T), MP.changecoefficienttype(AT, T)}
function Base.convert(::Type{BasicSemialgebraicSet{T, PT, AT}}, set::BasicSemialgebraicSet{T, PT, AT}) where {T, PT<:APL, AT<:AbstractAlgebraicSet}
    return set
end
function Base.convert(::Type{BasicSemialgebraicSet{T, PT, AT}}, set::BasicSemialgebraicSet) where {T, PT, AT}
    return BasicSemialgebraicSet{T, PT, AT}(set.V, set.p)
end

equalities(S::BasicSemialgebraicSet) = equalities(S.V)
addequality!(S::BasicSemialgebraicSet, p) = addequality!(S.V, p)
inequalities(S::BasicSemialgebraicSet) = S.p
addinequality!(S::BasicSemialgebraicSet, p) = push!(S.p, p)

Base.intersect(S::BasicSemialgebraicSet, T::BasicSemialgebraicSet) = BasicSemialgebraicSet(S.V ∩ T.V, [S.p; T.p])
Base.intersect(S::BasicSemialgebraicSet, T::AbstractAlgebraicSet) = BasicSemialgebraicSet(S.V ∩ T, copy(S.p))
Base.intersect(T::AbstractAlgebraicSet, S::BasicSemialgebraicSet) = intersect(S, T)
