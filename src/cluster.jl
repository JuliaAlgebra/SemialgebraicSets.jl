"""
    cluster_eigenvalues(_atol, v)

Clustering the values `v` following [CGT97].

[CGT97] Corless, R. M.; Gianni, P. M. & Trager, B. M.
*A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots*
Proceedings of the 1997 international symposium on Symbolic and algebraic computation, 1997, 133-140
"""
function cluster_eigenvalues(_atol, v)
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
            push!(λ, real(v[i]))
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

    # For eigenvalues not clustered yet, their eigenvalues is quite large.
    # Therefore, if we cluster all i, j close enough at once we might cluster too much
    # The technique used here is to cluster only the closest pair.
    # Once they are matched, a new atol is computed and if the cluster is complete,
    # this atol will be small which will avoid addition of new eigenvalues.
    while true
        σ = sortperm(λ)
        I = 0
        J = 0
        best = ONE
        for _i in eachindex(σ)
            i = σ[_i]
            for _j in 1:(_i-1)
                j = σ[_j]
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

    return clusters
end
