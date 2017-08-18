import Base: +

abstract type AbstractPolynomialIdeal end

struct PolynomialIdeal{T, PT<:APL{T}} <: AbstractPolynomialIdeal
    p::Vector{PT}
end
function PolynomialIdeal{T, PT}() where {T, PT<:APL{T}}
    PolynomialIdeal(PT[])
end
PolynomialIdeal(p::Vector{PT}) where {T, PT<:APL{T}} = PolynomialIdeal{T, PT}(p)
(+)(I::PolynomialIdeal, J::PolynomialIdeal) = PolynomialIdeal([I.p; J.p])
