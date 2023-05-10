export monomialbasis, ideal

abstract type AbstractPolynomialIdeal end

struct EmptyPolynomialIdeal <: AbstractPolynomialIdeal end
function ideal(p::FullSpace, algo = defaultgröbnerbasisalgorithm(p))
    return EmptyPolynomialIdeal()
end
Base.rem(p::AbstractPolynomialLike, I::EmptyPolynomialIdeal) = p

promote_for_division(::Type{T}) where {T} = T
promote_for_division(::Type{T}) where {T<:Integer} = Rational{big(T)}
promote_for_division(::Type{Rational{T}}) where {T} = Rational{big(T)}
function promote_for_division(::Type{Complex{T}}) where {T}
    return Complex{promote_for_division(T)}
end

mutable struct PolynomialIdeal{T,PT<:APL{T},A<:AbstractGröbnerBasisAlgorithm} <:
               AbstractPolynomialIdeal
    p::Vector{PT}
    gröbnerbasis::Bool
    algo::A
end
function PolynomialIdeal{T,PT}(
    p::Vector{PT},
    algo::A,
) where {T,PT<:APL{T},A<:AbstractGröbnerBasisAlgorithm}
    return PolynomialIdeal{T,PT,A}(p, false, algo)
end
function PolynomialIdeal(p::Vector{PT}, algo) where {T,PT<:APL{T}}
    S = promote_for_division(T)
    return PolynomialIdeal{S,polynomial_type(PT, S)}(
        AbstractVector{polynomial_type(PT, S)}(p),
        algo,
    )
end
function PolynomialIdeal{T,PT}() where {T,PT<:APL{T}}
    return PolynomialIdeal(PT[], defaultgröbnerbasisalgorithm(PT))
end

function Base.convert(
    ::Type{PolynomialIdeal{T,PT,A}},
    I::PolynomialIdeal{T,PT,A},
) where {T,PT<:APL{T},A<:AbstractGröbnerBasisAlgorithm}
    return I
end
function Base.convert(
    ::Type{PolynomialIdeal{T,PT,A}},
    I::PolynomialIdeal,
) where {T,PT,A}
    return PolynomialIdeal{T,PT,A}(I.p, I.gröbnerbasis, I.algo)
end

ideal(p, algo = defaultgröbnerbasisalgorithm(p)) = PolynomialIdeal(p, algo)

function Base.:+(I::PolynomialIdeal, J::PolynomialIdeal)
    return PolynomialIdeal([I.p; J.p], I.algo)
end

MP.variables(I::PolynomialIdeal) = variables(I.p)

function Base.rem(p::AbstractPolynomialLike, I::PolynomialIdeal)
    computegröbnerbasis!(I)
    return rem(p, I.p)
end

function computegröbnerbasis!(I::PolynomialIdeal)
    if !I.gröbnerbasis
        gröbnerbasis!(I.p, I.algo)
        I.gröbnerbasis = true
    end
end
function monomialbasis(I, vars = variables(I))
    computegröbnerbasis!(I)
    if isempty(I.p)
        return false, monomial_type(eltype(I.p))[]
    end
    mv = monomial_vector(leading_monomial.(I.p))
    # monomial_vector makes sure all monomials have the same variables
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
        return false, monomial_type(eltype(I.p))[]
    end
    return true,
    monomials(vars, 0:(sum(lv)-1), m -> !any(map(m2 -> divides(m2, m), mv)))
end
