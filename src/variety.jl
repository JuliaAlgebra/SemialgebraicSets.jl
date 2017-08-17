struct AlgebraicSet{T, PT<:APL{T}} <: AbstractAlgebraicSet
    I::PolynomialIdeal{T, PT}
end
function AlgebraicSet{T, PT}() where {T, PT<:APL{T}}
    AlgebraicSet(PolynomialIdeal{T, PT}())
end
AlgebraicSet(I::PolynomialIdeal{T, PT}) where {T, PT} = AlgebraicSet{T, PT}(I)
AlgebraicSet(p::Vector) = AlgebraicSet(PolynomialIdeal(p))

addequality!(V::AlgebraicSet, p) = push!(V.I.p, p)
Base.intersect(S::AlgebraicSet, T::AlgebraicSet) = AlgebraicSet(S.I + T.I)
