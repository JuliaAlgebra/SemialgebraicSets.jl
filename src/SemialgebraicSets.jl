module SemialgebraicSets

using Random

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
nequalities(::FullSpace) = 0
equalities(::FullSpace) = []
MP.changecoefficienttype(S::Type{FullSpace}, T::Type) = S

function MP.changecoefficienttype(set::AbstractSemialgebraicSet, T::Type)
    return convert(MP.changecoefficienttype(typeof(set), T), set)
end

Base.intersect(S::AbstractAlgebraicSet, T::FullSpace) = S
Base.intersect(S::FullSpace, T::AbstractAlgebraicSet) = T
Base.intersect(S::FullSpace, T::FullSpace) = S

# If `intersect(S, T)` is not implemented, this method will `StackOverflow`.
Base.intersect(S::AbstractSemialgebraicSet, T::AbstractSemialgebraicSet, args...; kws...) = intersect(intersect(S, T), args...; kws...)
# The keywords are only used when transforming `Element`
# into `BasicSemialgebraicSet`.
Base.intersect(set::AbstractSemialgebraicSet; kws...) = set

include("groebner.jl")
include("ideal.jl")
include("solve.jl")
include("variety.jl")
include("basic.jl")

include("fix.jl")
include("macro.jl")

end # module
