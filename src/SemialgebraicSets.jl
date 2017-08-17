module SemialgebraicSets

using MultivariatePolynomials

const APL = AbstractPolynomialLike

export AbstractSemialgebraicSet, AbstractBasicSemialgebraicSet, AbstractAlgebraicSet
export FullSpace, AlgebraicSet, BasicSemialgebraicSet, addequality!, addinequality!

# Semialgebraic set described by polynomials with coefficients in T
abstract type AbstractSemialgebraicSet end

abstract type AbstractBasicSemialgebraicSet <: AbstractSemialgebraicSet end
abstract type AbstractAlgebraicSet <: AbstractBasicSemialgebraicSet end

addinequality!(S::AbstractAlgebraicSet, p) = throw(ArgumentError("Cannot add inequality to an algebraic set"))

struct FullSpace <: AbstractAlgebraicSet
end

Base.intersect(S, T, U...) = intersect(intersect(S, T), U...)

include("ideal.jl")
include("variety.jl")
include("basic.jl")

include("macro.jl")

end # module
