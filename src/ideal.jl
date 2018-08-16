import Base: +

export monomialbasis, ideal

abstract type AbstractPolynomialIdeal end

struct EmptyPolynomialIdeal
end
ideal(p::FullSpace, algo=defaultgröbnerbasisalgorithm(p)) = EmptyPolynomialIdeal()
Base.rem(p::AbstractPolynomialLike, I::EmptyPolynomialIdeal) = p

mutable struct PolynomialIdeal{T, PT<:APL{T}, A<:AbstractGröbnerBasisAlgorithm} <: AbstractPolynomialIdeal
    p::Vector{PT}
    gröbnerbasis::Bool
    algo::A
end
function PolynomialIdeal{T, PT}(p::Vector{PT}, algo::A) where {T, PT<:APL{T}, A<:AbstractGröbnerBasisAlgorithm}
    PolynomialIdeal{T, PT, A}(p, false, algo)
end
function PolynomialIdeal(p::Vector{PT}, algo) where {T, PT<:APL{T}}
    S = Base.promote_op(/, T, T)
    PolynomialIdeal{S, polynomialtype(PT, S)}(AbstractVector{polynomialtype(PT, S)}(p), algo)
end
function PolynomialIdeal{T, PT}() where {T, PT<:APL{T}}
    PolynomialIdeal(PT[], defaultgröbnerbasisalgorithm(PT))
end

ideal(p, algo=defaultgröbnerbasisalgorithm(p)) = PolynomialIdeal(p, algo)

(+)(I::PolynomialIdeal, J::PolynomialIdeal) = PolynomialIdeal([I.p; J.p], I.algo)

MP.variables(I::PolynomialIdeal) = variables(I.p)

function Base.rem(p::AbstractPolynomialLike, I::PolynomialIdeal)
    computegröbnerbasis!(I)
    rem(p, I.p)
end

function computegröbnerbasis!(I::PolynomialIdeal)
    if !I.gröbnerbasis
        gröbnerbasis!(I.p, I.algo)
        I.gröbnerbasis = true
    end
end
function monomialbasis(I, vars=variables(I))
    computegröbnerbasis!(I)
    if isempty(I.p)
        return false, monomialtype(eltype(I.p))[]
    end
    mv = monovec(leadingmonomial.(I.p))
    # monovec makes sure all monomials have the same variables
    # if x_i^n is in the leading monomials, lv[i] <= n
    lv = -ones(Int, length(vars))
    for m in mv
        d = degree(m)
        if d == maximum(exponents(m))
            # univariate monomial
            if iszero(d)
                lv[:] .= 0
            else
                lv[argmax(exponents(m))] = d
            end
        end
    end
    if any(lv .< 0)
        return false, monomialtype(eltype(I.p))[]
    end
    return true, monomials(vars, 0:(sum(lv)-1), m -> !any(map(m2 -> divides(m2, m), mv)))
end
