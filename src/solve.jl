export algebraic_solver, ReorderedSchurMultiplicationMatricesSolver

"""
    AbstractAlgebraicSolver

Solver of algebraic equations.
"""
abstract type AbstractAlgebraicSolver end

function default_gröbner_basis_algorithm(p, ::AbstractAlgebraicSolver)
    return default_gröbner_basis_algorithm(p)
end

"""
    solvealgebraicequations(V::AbstractAlgebraicSet, algo::AbstractAlgebraicSolver)::Union{Nothing, Vector{<:Vector}}}

Solve the algebraic equations for which `V` is the set of solutions using the algorithm `algo`.
Returns a nullable which is `null` if `V` is not zero-dimensional and is the list of solutions otherwise.
"""
solve(V::AbstractAlgebraicSet) = solve(V, default_algebraic_solver(V))

"""
    AbstractMultiplicationMatricesAlgorithm

Algorithm computing multiplication matrices from algebraic equations.
"""
abstract type AbstractMultiplicationMatricesAlgorithm end

"""
    multiplication_matrices(V::AbstractAlgebraicSet, algo::AbstractMultiplicationMatricesAlgorithm)

Computing multiplication matrices from the algebraic equations for which `V` is the set of solution using the algorithm `algo`.
Returns a nullable which is `null` if `V` is not zero-dimensional and is the list of multiplication matrices otherwise.
"""
function multiplication_matrices end

"""
    AbstractMultiplicationMatricesSolver

Solver of algebraic equations using multiplication matrices.
"""
abstract type AbstractMultiplicationMatricesSolver end

struct SolverUsingMultiplicationMatrices{
    A<:AbstractMultiplicationMatricesAlgorithm,
    S<:AbstractMultiplicationMatricesSolver,
} <: AbstractAlgebraicSolver
    algo::A
    solver::S
end

function solve(V, solver::SolverUsingMultiplicationMatrices)
    Ms = multiplication_matrices(V, solver.algo)
    if Ms === nothing
        nothing
    else
        solve(Ms, solver.solver)
    end
end

struct MultiplicationMatrices{Ms}
    matrices::Ms
end

struct GröbnerBasisMultiplicationMatricesAlgorithm <:
       AbstractMultiplicationMatricesAlgorithm end

function multiplication_matrix(V, v::AbstractVariable, B)
    M = Matrix{eltype(eltype(V))}(undef, length(B), length(B))
    for i in 1:length(B)
        p = rem(v * B[i], equalities(V))
        M[:, i] = coefficients(p, B)
    end
    return M
end

function multiplication_matrices(
    V,
    ::GröbnerBasisMultiplicationMatricesAlgorithm,
)
    vars = variables(V.I)
    B = standard_monomials(V.I, vars)
    if isnothing(B)
        return
    else
        n = length(vars)
        if iszero(n)
            return Matrix{eltype(eltype(V))}[]
        else
            return MultiplicationMatrices([
                multiplication_matrix(V, v, B) for v in vars
            ])
        end
    end
end

include("schur.jl")

"""
Corless, R. M.; Gianni, P. M. & Trager, B. M. A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots Proceedings of the 1997 international symposium on Symbolic and algebraic computation, 1997, 133-140
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

function solve(
    Ms::MultiplicationMatrices,
    solver::ReorderedSchurMultiplicationMatricesSolver,
)
    λ = rand(solver.rng, length(Ms.matrices))
    λ /= sum(λ)
    return _solve_multiplication_matrices(Ms.matrices, λ, solver)
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

function algebraic_solver(
    algo::AbstractMultiplicationMatricesAlgorithm,
    solver::AbstractMultiplicationMatricesSolver,
)
    return SolverUsingMultiplicationMatrices(algo, solver)
end

function promote_for(
    ::Type{T},
    ::Type{<:SolverUsingMultiplicationMatrices},
) where {T}
    return float(T)
end

function default_multiplication_matrices_algorithm(p)
    return GröbnerBasisMultiplicationMatricesAlgorithm()
end
function default_multiplication_matrices_solver(::Type{T}) where {T}
    return ReorderedSchurMultiplicationMatricesSolver{T}()
end
function default_multiplication_matrices_solver(
    ::AbstractVector{PT},
) where {T,PT<:_APL{T}}
    return default_multiplication_matrices_solver(T)
end

function default_algebraic_solver(p)
    return algebraic_solver(
        default_multiplication_matrices_algorithm(p),
        default_multiplication_matrices_solver(p),
    )
end
function default_algebraic_solver(
    p,
    algo::AbstractMultiplicationMatricesAlgorithm,
)
    return algebraic_solver(algo, default_multiplication_matrices_solver(p))
end
function default_algebraic_solver(
    p,
    solver::AbstractMultiplicationMatricesSolver,
)
    return algebraic_solver(
        default_multiplication_matrices_algorithm(p),
        solver,
    )
end

function default_algebraic_solver(p, ::AbstractGröbnerBasisAlgorithm)
    return default_algebraic_solver(p)
end
