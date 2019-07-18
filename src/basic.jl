export inequalities, basicsemialgebraicset

struct BasicSemialgebraicSet{T, PT<:APL{T}, AT <: AbstractAlgebraicSet} <: AbstractBasicSemialgebraicSet
    V::AT
    p::Vector{PT}
end
function BasicSemialgebraicSet{T, PT}() where {T, PT<:APL{T}}
    BasicSemialgebraicSet(AlgebraicSet{T, PT}(), PT[])
end
function BasicSemialgebraicSet(V::AlgebraicSet{T, PT, A, S}, p::Vector{PT}) where {T, PT<:APL, A, S}
    BasicSemialgebraicSet{T, PT, typeof(V)}(V, p)
end
function BasicSemialgebraicSet(V::AlgebraicSet{T, PT, A, S}, p::Vector{PS}) where {T, PT<:APL, PS<:APL, A, S}
    PU = promote_type(PT, PS)
    BasicSemialgebraicSet(convert(AlgebraicSet{T, PU, A, S}, V), Vector{PU}(p))
end
#BasicSemialgebraicSet{T, PT<:APL{T}}(V::AlgebraicSet{T, PT}, p::Vector{PT}) = BasicSemialgebraicSet{T, PT}(V, p)
function basicsemialgebraicset(V, p)
    BasicSemialgebraicSet(V, p)
end

equalities(S::BasicSemialgebraicSet) = equalities(S.V)
addequality!(S::BasicSemialgebraicSet, p) = addequality!(S.V, p)
inequalities(S::BasicSemialgebraicSet) = S.p
addinequality!(S::BasicSemialgebraicSet, p) = push!(S.p, p)

Base.intersect(S::BasicSemialgebraicSet, T::BasicSemialgebraicSet) = intersect(BasicSemialgebraicSet(S.V ∩ T.V, [S.p; T.p]))
Base.intersect(S::BasicSemialgebraicSet, T::AbstractAlgebraicSet) = intersect(BasicSemialgebraicSet(S.V ∩ T, copy(S.p)))
Base.intersect(T::AbstractAlgebraicSet, S::BasicSemialgebraicSet) = intersect(S, T)
