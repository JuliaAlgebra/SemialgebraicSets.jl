module SemialgebraicSets

using Random

import MutableArithmetics
const MA = MutableArithmetics

using MultivariatePolynomials
const MP = MultivariatePolynomials

const APL = AbstractPolynomialLike

import CommonSolve: solve

export AbstractSemialgebraicSet,
    AbstractBasicSemialgebraicSet, AbstractAlgebraicSet
export FullSpace,
    AlgebraicSet, BasicSemialgebraicSet, add_equality!, add_inequality!

# Semialgebraic set described by polynomials with coefficients in T
abstract type AbstractSemialgebraicSet end

abstract type AbstractBasicSemialgebraicSet <: AbstractSemialgebraicSet end
abstract type AbstractAlgebraicSet <: AbstractBasicSemialgebraicSet end

function add_inequality!(S::AbstractAlgebraicSet, p)
    throw(ArgumentError("Cannot add inequality to an algebraic set"))
end

struct FullSpace <: AbstractAlgebraicSet end
function Base.show(io::IO, ::FullSpace)
    return print(io, "R^n")
end
nequalities(::FullSpace) = 0
equalities(::FullSpace) = []
MP.similar_type(S::Type{FullSpace}, T::Type) = S

function Base.similar(set::AbstractSemialgebraicSet, T::Type)
    return convert(MP.similar_type(typeof(set), T), set)
end

Base.intersect(S::AbstractAlgebraicSet, T::FullSpace) = S
Base.intersect(S::FullSpace, T::AbstractAlgebraicSet) = T
Base.intersect(S::FullSpace, T::FullSpace) = S

# If `intersect(S, T)` is not implemented, this method will `StackOverflow`.
function Base.intersect(
    S::AbstractSemialgebraicSet,
    T::AbstractSemialgebraicSet,
    args...;
    kws...,
)
    return intersect(intersect(S, T), args...; kws...)
end
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

include("deprecate.jl")

end # module
