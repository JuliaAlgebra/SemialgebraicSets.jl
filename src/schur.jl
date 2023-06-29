using LinearAlgebra

# If i, i+1 are conjugate pair, then they need to be either both in I or both not in I.
# If one of them is in I and the other is not then LAPACK will consider that both of them are in I.
function condition_number(sf::Schur, I)
    n = length(sf.values)
    select = zeros(LinearAlgebra.BlasInt, n)
    for i in I
        select[i] = 1
    end
    return LinearAlgebra.LAPACK.trsen!(
        'E',
        'N',
        select,
        copy(sf.T),
        copy(sf.Z),
    )[4]
end

# Manocha, D. & Demmel, J. Algorithms for intersecting parametric and algebraic curves II: multiple intersections
# Graphical Models and Image Processing, Elsevier, 1995, 57, 81-100
function clusterordschur(M::AbstractMatrix{<:Real}, ɛ)
    if isempty(M)
        # See bug JuliaLang/julia#...
        return Matrix{float(eltype(M))}(undef, 0, 0), Vector{Int}[]
    else
        return _clusterordschur(M, ɛ)
    end
end
function _clusterordschur(M::AbstractMatrix{BigFloat}, ɛ)
    Z, clusters = _clusterordschur(Float64.(M), ɛ)
    return convert(Matrix{BigFloat}, Z), clusters
end
function _clusterordschur(M::AbstractMatrix{<:Real}, ɛ)
    # M = Z * T * Z' and "values" gives the eigenvalues
    sf = LinearAlgebra.schur(M)
    Z = sf.Z
    v = sf.values
    # documentation says that the error on the eigenvalues is ɛ * norm(T) / condition_number
    nT = norm(sf.T)

    _atol(I) = ɛ * nT / condition_number(sf, I)

    A = typeof(_atol(1))
    V = real(eltype(v))

    ONE = one(one(V) / one(A))

    clusters = Vector{Int}[]
    λ = V[]
    atol = A[]
    # condition_number requires that conjugate pair need to be treated together so we first need to handle them
    # If they are in the same cluster then pair them, otherwise it is complex solution so we reject them
    i = firstindex(v)
    while i <= lastindex(v)
        if isreal(v[i])
            push!(clusters, [i])
            push!(λ, v[i])
            push!(atol, _atol(i))
            i += 1
        else
            @assert i < lastindex(v) && !isreal(v[i+1])
            pairatol = _atol([i, i + 1])
            if abs(v[i] - v[i+1]) / pairatol < ONE
                # Pair conjugate pairs into a real eigenvalue
                push!(clusters, [i, i + 1])
                push!(λ, real((v[i] + v[i+1]) / 2)) # The imaginary part should be zero anyway
                push!(atol, pairatol)
            end
            i += 2
        end
    end
    σ = sortperm(λ)

    # For eigenvalues not clustered yet, their eigenvalues is quite large.
    # Therefore, if we cluster all i, j close enough at once we might cluster too much
    # The technique used here is to cluster only the closest pair.
    # Once they are matched, a new atol is computed and if the cluster is complete,
    # this atol will be small which will avoid addition of new eigenvalues.
    while true
        I = 0
        J = 0
        best = ONE
        for i in eachindex(clusters)
            for j in 1:(i-1)
                d = abs(λ[i] - λ[j]) / min(atol[i], atol[j])
                if d < best
                    I = i
                    J = j
                    best = d
                end
            end
        end
        if best < ONE
            # merge I with J
            nI = length(clusters[I])
            nJ = length(clusters[J])
            λ[I] = (λ[I] * nI + λ[J] * nJ) / (nI + nJ)
            append!(clusters[I], clusters[J])
            atol[I] = _atol(clusters[I])
            deleteat!(λ, J)
            deleteat!(clusters, J)
            deleteat!(atol, J)
        else
            break
        end
    end
    return Z, clusters
end

"""
    struct ReorderedSchurMultiplicationMatricesSolver{T,RNGT<:Random.AbstractRNG} <:
        AbstractMultiplicationMatricesSolver
        ɛ::T
        rng::RNGT
    end

Simultaneous diagonalization of commuting matrices using the method of [CGT97].

[CGT97] Corless, R. M.; Gianni, P. M. & Trager, B. M.
*A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots*
Proceedings of the 1997 international symposium on Symbolic and algebraic computation, 1997, 133-140
"""
struct ReorderedSchurMultiplicationMatricesSolver{T,RNGT<:Random.AbstractRNG} <:
       AbstractMultiplicationMatricesSolver
    ɛ::T
    rng::RNGT
end
function ReorderedSchurMultiplicationMatricesSolver(ɛ)
    return ReorderedSchurMultiplicationMatricesSolver(ɛ, Random.GLOBAL_RNG)
end
function ReorderedSchurMultiplicationMatricesSolver{T}() where {T}
    return ReorderedSchurMultiplicationMatricesSolver(Base.rtoldefault(real(T)))
end

# Deterministic part
function _solve_multiplication_matrices(
    Ms::AbstractVector{<:AbstractMatrix{T}},
    λ,
    solver::ReorderedSchurMultiplicationMatricesSolver,
) where {T<:Real}
    @assert length(Ms) == length(λ)
    n = length(λ)
    Z, clusters = clusterordschur(sum(λ .* Ms), solver.ɛ)
    r = length(clusters)
    vals = [zeros(T, n) for k in 1:r]
    for k in 1:r
        nk = length(clusters[k])
        for j in clusters[k]
            q = Z[:, j]
            for i in 1:n
                vals[k][i] += dot(q, Ms[i] * q) / nk
            end
        end
    end
    return vals
end
