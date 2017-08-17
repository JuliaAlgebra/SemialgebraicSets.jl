struct BasicSemialgebraicSet{T, PT<:APL{T}} <: AbstractBasicSemialgebraicSet
    V::AlgebraicSet{T, PT}
    p::Vector{PT}
end
function BasicSemialgebraicSet{T, PT}() where {T, PT<:APL{T}}
    BasicSemialgebraicSet(AlgebraicSet{T, PT}(), PT[])
end
#BasicSemialgebraicSet{T, PT<:APL{T}}(V::AlgebraicSet{T, PT}, p::Vector{PT}) = BasicSemialgebraicSet{T, PT}(V, p)

addequality!(S::BasicSemialgebraicSet, p) = addequality!(S.V, p)
addinequality!(S::BasicSemialgebraicSet, p) = push!(S.p, p)

Base.intersect(S::BasicSemialgebraicSet, T::BasicSemialgebraicSet) = BasicSemialgebraicSet(S.V ∩ T.V, [S.p; T.p])
Base.intersect(S::BasicSemialgebraicSet, T::AlgebraicSet) = BasicSemialgebraicSet(S.V ∩ T, copy(S.p))
Base.intersect(T::AlgebraicSet, S::BasicSemialgebraicSet) = intersect(S, T)
