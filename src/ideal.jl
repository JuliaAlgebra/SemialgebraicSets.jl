export standard_monomials, ideal

abstract type AbstractPolynomialIdeal end

struct EmptyPolynomialIdeal <: AbstractPolynomialIdeal end
function ideal(p::FullSpace, _ = default_gröbner_basis_algorithm(p))
    return EmptyPolynomialIdeal()
end
Base.rem(p::MP.AbstractPolynomialLike, ::EmptyPolynomialIdeal) = p

mutable struct PolynomialIdeal{
    T,
    PT<:_APL{T},
    A<:AbstractGröbnerBasisAlgorithm,
} <: AbstractPolynomialIdeal
    p::Vector{PT}
    gröbner_basis::Bool
    algo::A
end
MP.coefficient_type(::Type{<:PolynomialIdeal{T}}) where {T} = T

function PolynomialIdeal{T,PT}(
    p::Vector{PT},
    algo::A,
) where {T,PT<:_APL{T},A<:AbstractGröbnerBasisAlgorithm}
    return PolynomialIdeal{T,PT,A}(p, false, algo)
end
function PolynomialIdeal(p::Vector{PT}, algo) where {PT<:_APL}
    S = promote_for(MP.coefficient_type(PT), typeof(algo))
    return PolynomialIdeal{S,MP.polynomial_type(PT, S)}(
        AbstractVector{MP.polynomial_type(PT, S)}(p),
        algo,
    )
end
function PolynomialIdeal{T,PT}() where {T,PT<:_APL{T}}
    return PolynomialIdeal(PT[], default_gröbner_basis_algorithm(PT))
end

function Base.convert(
    ::Type{PolynomialIdeal{T,PT,A}},
    I::PolynomialIdeal{T,PT,A},
) where {T,PT<:_APL{T},A<:AbstractGröbnerBasisAlgorithm}
    return I
end
function Base.convert(
    ::Type{PolynomialIdeal{T,PT,A}},
    I::PolynomialIdeal,
) where {T,PT,A}
    return PolynomialIdeal{T,PT,A}(I.p, I.gröbner_basis, I.algo)
end

ideal(p, algo = default_gröbner_basis_algorithm(p)) = PolynomialIdeal(p, algo)

function Base.:+(I::PolynomialIdeal, J::PolynomialIdeal)
    return PolynomialIdeal(
        _promote_polys_to_common_variables([I.p; J.p]),
        I.algo,
    )
end

MP.variables(I::PolynomialIdeal) = MP.variables(I.p)

function Base.rem(p::MP.AbstractPolynomialLike, I::PolynomialIdeal)
    compute_gröbner_basis!(I)
    return rem(p, I.p)
end

function compute_gröbner_basis!(I::PolynomialIdeal)
    if !I.gröbner_basis
        gröbner_basis!(I.p, I.algo)
        I.gröbner_basis = true
    end
end

"""
    is_gröbner_basis(I::AbstractPolynomialIdeal)

Return whether the generators of `I` are known to form a Gröbner basis.
If it returns `true` then `rem(p, I)` does not need to compute a Gröbner basis
(it may still normalize the generators, e.g., make them monic).
If it returns `false`, `rem(p, I)` would first compute a Gröbner basis with
[`compute_gröbner_basis!`](@ref), which may be expensive.
"""
function is_gröbner_basis end

is_gröbner_basis(::EmptyPolynomialIdeal) = true

function is_gröbner_basis(I::PolynomialIdeal)
    # A single polynomial always forms a Gröbner basis
    return I.gröbner_basis || length(I.p) <= 1
end

"""
    standard_monomials(I::AbstractPolynomialIdeal, vars = MP.variables(I))

Return the vector of *standard monomials* in the variables `var` of the ideal
`I`. A monomial is a *standard monomial* if it is not the leading monomial of
any polynomial of the ideal. If there are infinitely such monomials, `nothing`
is returned.

See [CLO5, p. 38], where it is also called *basis monomial*, for more
information.

[CLO05] Cox, A. David & Little, John & O'Shea, Donal
*Using Algebraic Geometry*.
Graduate Texts in Mathematics, **2005**.
https://doi.org/10.1007/b138611
"""
function standard_monomials end

"""
    standard_monomials(I::PolynomialIdeal, vars, maxdegree::Integer)

Return the vector of *standard monomials* of the ideal `I` in the variables
`vars` of degree at most `maxdegree`. As the standard monomials form an order
ideal, the standard monomials of degree `d + 1` are obtained by multiplying
the standard monomials of degree `d` by each variable and keeping the standard
ones. The time complexity is therefore proportional to the number of standard
monomials and not to the number of monomials of degree at most `maxdegree`.
"""
function standard_monomials(I::PolynomialIdeal, vars, maxdegree::Integer)
    compute_gröbner_basis!(I)
    leading = [MP.exponents(MP.leading_monomial(p), vars) for p in I.p]
    is_standard(e) = !any(l -> all(l .<= e), leading)
    current = [zeros(Int, length(vars))]
    if maxdegree < 0 || !is_standard(first(current))
        empty!(current)
    end
    exponents = copy(current)
    for _ in 1:maxdegree
        next = Set{Vector{Int}}()
        for e in current
            for i in eachindex(vars)
                f = copy(e)
                f[i] += 1
                if is_standard(f)
                    push!(next, f)
                end
            end
        end
        current = collect(next)
        append!(exponents, current)
    end
    M = MP.monomial_type(prod(vars))
    return MP.monomial_vector(M[prod(vars .^ e) for e in exponents])
end

function standard_monomials(I, vars = MP.variables(I))
    compute_gröbner_basis!(I)
    if isempty(I.p)
        return nothing # Nonzero dimensional
    end
    mv = MP.monomial_vector(MP.leading_monomial.(I.p))
    # monomial_vector makes sure all monomials have the same variables
    # if x_i^n is in the leading monomials, lv[i] <= n
    lv = -ones(Int, length(vars))
    for m in mv
        d = MP.degree(m)
        if d == maximum(MP.exponents(m))
            # univariate monomial
            if iszero(d)
                lv[:] .= 0
            else
                lv[argmax(MP.exponents(m))] = d
            end
        end
    end
    if any(Base.Fix2(<, 0), lv)
        return nothing # Nonzero dimensional
    end
    return MP.monomials(
        vars,
        0:(sum(lv)-1),
        m -> !any(Base.Fix2(MP.divides, m), mv),
    )
end
