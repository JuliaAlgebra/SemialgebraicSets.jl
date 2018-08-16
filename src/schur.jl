import Compat.LinearAlgebra: Schur, BlasInt, checksquare, chkstride1
import Compat.LinearAlgebra.LAPACK: liblapack, chklapackerror, @blasfunc

using Compat.LinearAlgebra

if VERSION <= v"0.7-"
    # Taken from JuliaLang/julia/base/linalg/lapack.jl and fixed for compq == 'E'
    # See https://github.com/JuliaLang/julia/pull/23521
    for (trexc, trsen, tgsen, elty) in
        ((:dtrexc_, :dtrsen_, :dtgsen_, :Float64),
         (:strexc_, :strsen_, :stgsen_, :Float32))
        @eval begin
            # *     .. Scalar Arguments ..
            #       CHARACTER          COMPQ, JOB
            #       INTEGER            INFO, LDQ, LDT, LIWORK, LWORK, M, N
            #       DOUBLE PRECISION   S, SEP
            # *     ..
            # *     .. Array Arguments ..
            #       LOGICAL            SELECT( * )
            #       INTEGER            IWORK( * )
            #       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ), WR( * )
            function _trsen!(compq::Char, job::Char, select::StridedVector{BlasInt},
                            T::StridedMatrix{$elty}, Q::StridedMatrix{$elty})
                chkstride1(T, Q, select)
                n = checksquare(T)
                ldt = max(1, stride(T, 2))
                ldq = max(1, stride(Q, 2))
                wr = similar(T, $elty, n)
                wi = similar(T, $elty, n)
                m = sum(select)
                work = Vector{$elty}(1)
                lwork = BlasInt(-1)
                iwork = Vector{BlasInt}(1)
                liwork = BlasInt(-1)
                info = Ref{BlasInt}()
                select = convert(Array{BlasInt}, select)
                s = Ref{$elty}(zero($elty))
                sep = Ref{$elty}(zero($elty))
                for i = 1:2
                    ccall((@blasfunc($trsen), liblapack), Void,
                        (Ptr{UInt8}, Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt},
                        Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
                        Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt}, Ref{$elty}, Ref{$elty},
                        Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                        Ptr{BlasInt}),
                        &compq, &job, select, &n,
                        T, &ldt, Q, &ldq,
                        wr, wi, &m, s, sep,
                        work, &lwork, iwork, &liwork,
                        info)
                    chklapackerror(info[])
                    if i == 1 # only estimated optimal lwork, liwork
                        lwork  = BlasInt(real(work[1]))
                        liwork = BlasInt(real(iwork[1]))
                        work   = Vector{$elty}(lwork)
                        iwork  = Vector{BlasInt}(liwork)
                    end
                end
                T, Q, iszero(wi) ? wr : complex.(wr, wi), s[], sep[]
            end
        end
    end
else
    _trsen! = LinearAlgebra.LAPACK.trsen!
end

# If i, i+1 are conjugate pair, then they need to be either both in I or both not in I.
# If one of them is in I and the other is not then LAPACK will consider that both of them are in I.
function conditionnumber(sf::Schur, I)
    n = length(sf.values)
    select = zeros(BlasInt, n)
    for i in I
        select[i] = 1
    end
    _trsen!('E', 'N', select, copy(sf.T), copy(sf.Z))[4]
end

# Manocha, D. & Demmel, J. Algorithms for intersecting parametric and algebraic curves II: multiple intersections
# Graphical Models and Image Processing, Elsevier, 1995, 57, 81-100
function clusterordschur(M::AbstractMatrix{<:Real}, ɛ)
    if isempty(M)
        # See bug JuliaLang/julia#...
        return Matrix{float(eltype(M))}(undef, 0, 0), Vector{Int}[]
    else
        _clusterordschur(M, ɛ)
    end
end
function _clusterordschur(M::AbstractMatrix{<:Real}, ɛ)
    # M = Z * T * Z' and "values" gives the eigenvalues
    if VERSION >= v"0.7-"
        sf = schur(M)
        Z = sf.Z
        v = sf.values
    else
        sf = schurfact(M)
        Z = sf[:Z]
        v = sf[:values]
    end
    # documentation says that the error on the eigenvalues is ɛ * norm(T) / conditionnumber
    nT = norm(sf.T)

    _atol(I) = ɛ * nT / conditionnumber(sf, I)

    A = typeof(_atol(1))
    V = real(eltype(v))

    ONE = one(one(V) / one(A))

    clusters = Vector{Int}[]
    λ = V[]
    atol = A[]
    # conditionnumber requires that conjugate pair need to be treated together so we first need to handle them
    # If they are in the same cluster then pair them, otherwise it is complex solution so we reject them
    i = 1
    while i <= lastindex(v)
        if isreal(v[i])
            push!(clusters, [i])
            push!(λ, v[i])
            push!(atol, _atol(i))
            i += 1
        else
            @assert i < lastindex(v) && !isreal(v[i+1])
            pairatol = _atol([i, i+1])
            if abs(v[i] - v[i+1]) / pairatol < ONE
                # Pair conjugate pairs into a real eigenvalue
                push!(clusters, [i, i+1])
                push!(λ, real((v[i] + v[i+1]) / 2)) # The imaginary part should be zero anyway
                push!(atol, pairatol)
            end
            i += 2
        end
    end
    σ = sortperm(λ)

    # For eigenvalues not clustered yet, their eigenvalues is quite large.
    # Therefore, if we cluster all i, j close enough at once we migth cluster too much
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
    Z, clusters
end
