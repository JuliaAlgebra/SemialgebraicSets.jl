"""
    spolynomial(p::AbstractPolynomialLike, q::AbstractPolynomialLike)

Computes the S-polynomial of `p` and `q` defined by
```math
S(p, q) =  \\frac{m}{\\mathrm{\\mathsc{LT}}(p)} \\cdot p - \\frac{m}{\\mathrm{\\mathsc{LT}}(q)} \\cdot q
```
where ``m = \\mathrm{lcm}(\\mathrm{\\mathsc{LM}}(p), \\mathrm{\\mathsc{LM}}(q))``.
"""
function spolynomial(p::APL, q::APL)
    m = lcm(leadingmonomial(p), leadingmonomial(q))
    MultivariatePolynomials._div(m, leadingterm(p)) * p - MultivariatePolynomials._div(m, leadingterm(q)) * q
end

# Taken from Theorem 9 of
# Ideals, Varieties, and Algorithms
# Cox, Little and O'Shea, Fourth edition
function gröbnerbasis!(F::Vector{<:APL})
    n = length(F)
    B = Set{NTuple{2, Int}}(Iterators.filter(t -> t[1] < t[2], (i, j) for i in 1:n, j in 1:n))
    while !isempty(B)
        (i, j) = first(B)
        lmi = leadingmonomial(F[i])
        lmj = leadingmonomial(F[j])
        if !isconstant(gcd(lmi, lmj))
        end
    end
end
gröbnerbasis(F::Vector{<:APL}) = gröbnerbasis!(copy(F))
const groebnerbasis = gröbnerbasis
