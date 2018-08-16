__precompile__()

module SemialgebraicSets

using Compat
using Compat.Random

using MultivariatePolynomials
const MP = MultivariatePolynomials

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
function Base.show(io::IO, ::FullSpace)
    print(io, "R^n")
end

Base.intersect(S::AbstractSemialgebraicSet, T::AbstractSemialgebraicSet, U::AbstractSemialgebraicSet...) = intersect(intersect(S, T), U...)

include("groebner.jl")
include("ideal.jl")
include("solve.jl")
include("variety.jl")
include("basic.jl")

include("macro.jl")

end # module
