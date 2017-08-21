import Base: +

abstract type AbstractPolynomialIdeal end

mutable struct PolynomialIdeal{T, PT<:APL{T}} <: AbstractPolynomialIdeal
    p::Vector{PT}
    gröbnerbasis::Bool
end
function PolynomialIdeal{T, PT}(p::Vector{PT}) where {T, PT<:APL{T}}
    PolynomialIdeal{T, PT}(p, false)
end
function PolynomialIdeal(p::Vector{PT}) where {T, PT<:APL{T}}
    S = Base.promote_op(/, T, T)
    PolynomialIdeal{S, polynomialtype(PT, S)}(AbstractVector{polynomialtype(PT, S)}(p))
end
function PolynomialIdeal{T, PT}() where {T, PT<:APL{T}}
    PolynomialIdeal(PT[])
end
(+)(I::PolynomialIdeal, J::PolynomialIdeal) = PolynomialIdeal([I.p; J.p])

function computegröbnerbasis!(I::PolynomialIdeal)
    if !I.gröbnerbasis
        gröbnerbasis!(I.p)
        I.gröbnerbasis = true
    end
end
function monomialbasis(I)
    computegröbnerbasis!(I)
    if isempty(I.p)
        return false, monomialtype(eltype(I.p))[]
    end
    mv = monovec(leadingmonomial.(I.p))
    # monovec makes sure all monomials have the same variables
    # if x_i^n is in the leading monomials, lv[i] <= n
    lv = -ones(Int, nvariables(mv))
    for m in mv
        d = deg(m)
        if d == maximum(exponents(m))
            # univariate monomial
            if iszero(d)
                lv[:] = 0
            else
                lv[indmax(exponents(m))] = d
            end
        end
    end
    if any(lv .< 0)
        return false, monomialtype(eltype(I.p))[]
    end
    return true, monomials(variables(I.p), 0:(sum(lv)-1), m -> !any(map(m2 -> divides(m2, m), mv)))
end
