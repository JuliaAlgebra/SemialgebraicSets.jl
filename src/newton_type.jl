# This file is largely inspired from Bernard Mourrain's MultivariateSeries/diagonalization.jl

export NewtonTypeDiagonalization

# norm of off diagonal terms of a square matrix
function norm_off(M)
    n = LinearAlgebra.checksquare(M)
    if n > 1
        return sqrt(
            sum(abs2(M[i, j]) + abs2(M[j, i]) for i in 1:n for j in i+1:n),
        )
    else
        return zero(eltype(M))
    end
end

"""
    diagonalization_iter(D::Vector{<:AbstractMatrix{T}}) where {T}

Given the vector `D` of `[F_i * E_i, F_i * M_1 * E_i, ..., F_i * M_p * E_i]`,
computes the matrix `X` (resp. `Y`) corresponding to ``(I_n + X_i)``
(resp. ``(I_n + Y_i)``) of [KMY22, Theorem 5] that solves equations
[KMY22, (26)-(29)].
"""
function diagonalization_iter(D::Vector{<:AbstractMatrix{T}}) where {T}
    n = LinearAlgebra.checksquare(D[1])
    s = length(D)

    X = fill(zero(T), n, n)
    Y = fill(zero(T), n, n)

    A = fill(zero(T), s, 2)
    b = fill(zero(T), s)
    for i in 1:n
        for j in 1:n
            if i != j
                for k in 1:s
                    A[k, 1] = D[k][i, i]
                    A[k, 2] = D[k][j, j]
                    b[k] = -D[k][i, j]
                end
                v = A \ b
                X[i, j] = v[1]
                Y[i, j] = v[2]
            end
        end
    end
    for i in 1:n
        X[i, i] = 1
        Y[i, i] = 1
    end
    return X, Y
end

"""
    struct NewtonTypeDiagonalization{T,RNGT} <: AbstractMultiplicationMatricesSolver
        max_iter::Int
        ε::T
        tol::T
        rng::RNGT
    end

Simultaneous diagonalization of commuting matrices using the method of [KMY22, Theorem 5].

[KMY22] Khouja, Rima, Mourrain, Bernard, and Yakoubsohn, Jean-Claude.
*Newton-type methods for simultaneous matrix diagonalization.*
Calcolo 59.4 (2022): 38.
"""
struct NewtonTypeDiagonalization{T,RNGT} <: AbstractMultiplicationMatricesSolver
    max_iter::Int
    ε::T
    tol::T
    rng::RNGT
end
# These were the values in MultivariateSeries/diagonalization.jl
function NewtonTypeDiagonalization(max_iter, ε, tol)
    return NewtonTypeDiagonalization(max_iter, ε, tol, Random.GLOBAL_RNG)
end
NewtonTypeDiagonalization() = NewtonTypeDiagonalization(10, 1e-3, 5e-2)

function _eigvecs(M::AbstractMatrix{BigFloat})
    ev = _eigvecs(Float64.(M))
    return convert(Matrix{BigFloat}, ev)
end
# `eigvecs` is failing some tests with a non-invertible `ev`
#_eigvecs(M::AbstractMatrix) = LinearAlgebra.eigvecs(M)
_eigvecs(M::AbstractMatrix) = LinearAlgebra.schur(M).vectors

function _solve_multiplication_matrices(M, λ, solver::NewtonTypeDiagonalization)
    @assert length(M) == length(λ)
    n = length(λ)
    r = LinearAlgebra.checksquare(M[1])

    M1 = sum(λ .* M)
    E = _eigvecs(M1)

    # With `eigvecs`, we should do `inv` but with `schur` we can just transpose
    #F = inv(E)
    F = E'

    D = vcat(
        # Add one matrix for the equation `F_i * E_i = I` 
        # constraining `E_i` to be invertible
        [Matrix{eltype(M[1])}(I, r, r)],
        [F * M[i] * E for i in eachindex(M)],
    )
    err = sum(norm_off.(D))
    Δ = sum(norm.(D))

    nit = 0

    if err / Δ > solver.tol
        Δ = err
        while nit < solver.max_iter && Δ > solver.ε
            err0 = err
            X, Y = diagonalization_iter(D)
            # From [KMY22, Theorem 5]
            # Z_{i,k} + ∑_{i,k}
            # = F_i * M_k * E_i
            # = (I_n * Y_i) * (F_{i-1} * M_k * E_{i-1}) * (I_n * X_i)
            # = (I_n * Y_i) * D[i] * (I_n * X_i)
            # = Y * D[i] * X
            D = [Y * D[i] * X for i in eachindex(D)]
            # E_{i+1} = E_i * (I_n * X_i) from [KMY22, Theorem 5]
            E = E * X
            # F_{i+1} = (I_n * Y_i) * F_i from [KMY22, Theorem 5]
            F = Y * F
            nit += 1
            err = sum(norm_off.(D))
            Δ = err0 - err
        end
    end

    return [[D[j+1][i, i] / D[1][i, i] for j in 1:n] for i in 1:r]

    #    # I implemented this when I was analysing the result after only zero iteration so unsure if it's useful
    #    d = LinearAlgebra.Diagonal(sqrt.(inv.(LinearAlgebra.diag(D[1]))))
    #    Λ = [d * D[j+1] * d for j in 1:n]

    #    # `Λ` can be decomposed into blocks corresponding to the same eigenvalue
    #    # These blocks can have off-diagonal entries so we further diagonalize
    #    # with
    #    # FIXME `solver.tol` or `solver.ε` or ?
    #    sub_solver = ReorderedSchurMultiplicationMatricesSolver(solver.ε, solver.rng)
    #    sols = Vector{eltype(Λ[1])}[]
    #    i = 1
    #    while i <= r
    #        j = findfirst((i+1):r) do j
    #            all(1:n) do k
    #                # FIXME `solver.tol` or `solver.ε` or ?
    #                return !isapprox(Λ[k][j, j], Λ[k][i, i], rtol = solver.ε)
    #            end
    #        end
    #        j = something(j, r - i + 1)
    #        I = i:(i+j-1)
    #        sub_matrices = MultiplicationMatrices([Λ[j][I, I] for j in 1:n])
    #        append!(sols, solve(sub_matrices, sub_solver))
    #        i += j
    #    end
    #    return sols
end
