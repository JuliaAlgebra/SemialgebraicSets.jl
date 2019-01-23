export inequalities, basicsemialgebraicset

struct BasicSemialgebraicSet{T, PT<:APL{T}, AT <: AlgebraicSet{T, PT}} <: AbstractBasicSemialgebraicSet
    V::AT
    p::Vector{PT}
end
function BasicSemialgebraicSet{T, PT}() where {T, PT<:APL{T}}
    BasicSemialgebraicSet(AlgebraicSet{T, PT}(), PT[])
end
function BasicSemialgebraicSet(V::AlgebraicSet{T, PT, A, S}, p::Vector{PS}) where {T, PT, PS, A, S}
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

Base.intersect(S::BasicSemialgebraicSet, T::BasicSemialgebraicSet) = BasicSemialgebraicSet(S.V ∩ T.V, [S.p; T.p])
Base.intersect(S::BasicSemialgebraicSet, T::AlgebraicSet) = BasicSemialgebraicSet(S.V ∩ T, copy(S.p))
Base.intersect(T::AlgebraicSet, S::BasicSemialgebraicSet) = intersect(S, T)
