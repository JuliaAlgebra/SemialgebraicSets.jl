export spolynomial, gröbnerbasis, groebnerbasis
export presort!, dummyselection, normalselection

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

function reducebasis!(F::AbstractVector{<:APL})
    changed = true
    while changed
        changed = false
        keep = collect(eachindex(F))
        del = Int[]
        for j in eachindex(F)
            G = @view F[setdiff(keep, j)]
            r = rem(F[j], G)
            if isapproxzero(r)
                deleteat!(keep, findfirst(j))
                push!(del, j)
                changed = true
            elseif leadingmonomial(r) != leadingmonomial(F[j])
                F[j] = r
                changed = true
            end
        end
        deleteat!(F, del)
    end
end

clean!(F::AbstractVector{<:APL}) = map!(monic, F, F)

ext(i, j) = (min(i, j), max(i, j))

lcmlm(F, i, j) = lcm(leadingmonomial(F[i]), leadingmonomial(F[j]))

function criterion(F, B, i, j)
    m = lcmlm(F, i, j)
    for k in eachindex(F)
        if k != i && k != j &&
            !(ext(i, k) in B) && !(ext(j, k) in B) &&
            divides(leadingmonomial(F[k]), m)
            return true
        end
    end
    return false
end

# Preprocessing of F

# "In practice, some effort could be saved on average if we arranged the f_i so that their
# leading terms are listed in increasing order with respect to the chosen monomial ordering.
# Since the smaller LT(fi) are more likely to be used during the division algorithm,
# listing them earlier means that fewer comparisons will be required."
#
# Ideals, Varieties, and Algorithms
# Cox, Little and O'Shea, Fourth edition p. 116
function presort!(F)
    sort!(F, by=leadingmonomial)
end

# Selection of (i, j) ∈ B

function dummyselection(F, B)
    first(B)
end
# "BUCHBERGER (1985) suggests that there will often be some savings if we pick
# (i, j) ∈ B such that lcm(LM(fi), LM(fj)) is as small as possible.
# The corresponding S-polynomials will tend to yield any nonzero remainders
# (and new elements of the Gröbner basis) sooner in the process,
# so there will be more of a chance that subsequent remainders rem(S(fi, fj), G) will be zero.
# This normal selection strategy is discussed in more detail in BECKER and WEISPFENNING (1993),
# BUCHBERGER (1985) and GEBAUER and MÖLLER (1988)."
#
# Ideals, Varieties, and Algorithms
# Cox, Little and O'Shea, Fourth edition p. 116
function normalselection(F, B)
    best = first(B)
    m = lcmlm(F, best...)
    for (i, j) in Iterators.drop(B, 1)
        cur = lcmlm(F, i, j)
        if cur < m
            best = (i, j)
            m = cur
        end
    end
    best
end
# TODO sugar and double sugar selections

struct Buchberger
    pre!::Function
    sel::Function
end
#Buchberger() = Buchberger(presort!, dummyselection)
Buchberger() = Buchberger(presort!, normalselection)

# Taken from Theorem 9 of
# Ideals, Varieties, and Algorithms
# Cox, Little and O'Shea, Fourth edition
function gröbnerbasis!(F::AbstractVector{<:APL}, algo=Buchberger())
    algo.pre!(F)
    B = Set{NTuple{2, Int}}(Iterators.filter(t -> t[1] < t[2], (i, j) for i in eachindex(F), j in eachindex(F)))
    count = 0
    while !isempty(B)
        (i, j) = algo.sel(F, B)
        lmi = leadingmonomial(F[i])
        lmj = leadingmonomial(F[j])
        if !isconstant(gcd(lmi, lmj)) && !criterion(F, B, i, j)
            count += 1
            r = rem(spolynomial(F[i], F[j]), F)
            if !isapproxzero(r)
                I = eachindex(F)
                push!(F, r)
                n = last(eachindex(F))
                for i in I
                    push!(B, (i, n))
                end
            end
        end
        delete!(B, (i, j))
    end
    @show count
    reducebasis!(F)
    clean!(F)
    F
end
function gröbnerbasis(F::Vector{<:APL})
    T = Base.promote_op(rem, eltype(F), typeof(F))
    gröbnerbasis!(copy!(similar(F, T), F))
end
const groebnerbasis = gröbnerbasis
