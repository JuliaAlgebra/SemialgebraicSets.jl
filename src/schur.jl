import Base.LinAlg: Schur, BlasInt, checksquare, chkstride1
import Base.LAPACK: liblapack, chklapackerror, @blasfunc

# Taken from JuliaLang/julia/base/linalg/lapack.jl
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
function conditionnumber(sf::Schur, i)
    n = length(sf.values)
    select = zeros(BlasInt, n)
    select[i] = 1
    _trsen!('E', 'N', select, copy(sf.T), copy(sf.Z))[4]
end

# Algorithms for intersecting parametric and algebraic curves II: multiple intersections
# Manocha, Dinesh and Demmel, James
# Graphical Models and Image Processing, 1995
function clusterordschur(M, ɛ)
    sf = schurfact(M)
    # M = Z * T * Z' and "values" gives the eigenvalues
    Z = sf[:Z]
    v = copy(sf[:values])
    # documentation says that the error on the eigenvalues is ɛ * norm(T) / conditionnumber
    nT = norm(sf.T)
    _atol(i) = ɛ * nT / conditionnumber(sf, i)
    n = length(v)
    atol = _atol.(1:n)
    # Clustering
    clusters = Vector{Int}[[i] for i in 1:n]
    ONE = abs(one(Base.promote_op(/, eltype(v), eltype(atol))))
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
                d = abs(v[i] - v[j]) / min(atol[i], atol[j])
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
            v[I] = (v[I] * nI + v[J] * nJ) / (nI + nJ)
            append!(clusters[I], clusters[J])
            atol[I] = _atol(clusters[I])
            deleteat!(v, J)
            deleteat!(clusters, J)
            deleteat!(atol, J)
        else
            break
        end
    end
    Z, clusters[map(k -> isapproxzero(imag(v[k]); ztol=atol[k]), eachindex(clusters))]
end
