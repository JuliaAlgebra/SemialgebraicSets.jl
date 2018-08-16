export algebraicsolver, ReorderedSchurMultiplicationMatricesSolver

"""
    AbstractAlgebraicSolver

Solver of algebraic equations.
"""
abstract type AbstractAlgebraicSolver end

"""
    solvealgebraicequations(V::AbstractAlgebraicSet, algo::AbstractAlgebraicSolver)::Union{Nothing, Vector{<:Vector}}}

Solve the algebraic equations for which `V` is the set of solutions using the algorithm `algo`.
Returns a nullable which is `null` if `V` is not zero-dimensional and is the list of solutions otherwise.
"""
solvealgebraicequations(V::AbstractAlgebraicSet) = solvealgebraicequations(V, defaultalgebraicsolver(V))

"""
    AbstractMultiplicationMatricesAlgorithm

Algorithm computing multiplication matrices from algebraic equations.
"""
abstract type AbstractMultiplicationMatricesAlgorithm end

"""
    multiplicationmatrices(V::AbstractAlgebraicSet, algo::AbstractMultiplicationMatricesAlgorithm)::Union{Nothing, Vector{<:AbstractMatrix}}

Computing multiplication matrices from the algebraic equations for which `V` is the set of solution using the algorithm `algo`.
Returns a nullable which is `null` if `V` is not zero-dimensional and is the list of multiplication matrices otherwise.
"""
function multiplicationmatrices end

"""
    AbstractMultiplicationMatricesSolver

Solver of algebraic equations using multiplication matrices.
"""
abstract type AbstractMultiplicationMatricesSolver end

"""
    solvemultiplicationmatrices(Ms::AbstractVector{<:AbstractMatrix{T}}, algo::AbstractMultiplicationMatricesSolver)::Vector{Vector{T}} where T

Solve the algebraic equations having multiplication matrices `Ms` using the algorithm `algo`.
Returns the list of solutions.
"""
function solvemultiplicationmatrices end


struct SolverUsingMultiplicationMatrices{A<:AbstractMultiplicationMatricesAlgorithm, S<:AbstractMultiplicationMatricesSolver} <: AbstractAlgebraicSolver
    algo::A
    solver::S
end

function solvealgebraicequations(V::AbstractAlgebraicSet, solver::SolverUsingMultiplicationMatrices)
    Ms = multiplicationmatrices(V, solver.algo)
    if Ms === nothing
        nothing
    else
        solvemultiplicationmatrices(Ms, solver.solver)
    end
end

struct GröbnerBasisMultiplicationMatricesAlgorithm <: AbstractMultiplicationMatricesAlgorithm
end

function multiplicationmatrix(V::AbstractAlgebraicSet, v::AbstractVariable, B)
    M = Matrix{eltype(eltype(V))}(undef, length(B), length(B))
    for i in 1:length(B)
        p = rem(v * B[i], equalities(V))
        M[:, i] = coefficients(p, B)
    end
    M
end

function multiplicationmatrices(V::AbstractAlgebraicSet, algo::GröbnerBasisMultiplicationMatricesAlgorithm)
    vars = variables(V.I)
    iszd, B = monomialbasis(V.I, vars)
    if !iszd
        nothing
    else
        n = length(vars)
        if iszero(n)
            Matrix{eltype(eltype(T))}[]
        else
            [multiplicationmatrix(V, v, B) for v in vars]
        end
    end
end

include("schur.jl")

"""
Corless, R. M.; Gianni, P. M. & Trager, B. M. A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots Proceedings of the 1997 international symposium on Symbolic and algebraic computation, 1997, 133-140
"""
struct ReorderedSchurMultiplicationMatricesSolver{T, RNGT <: Compat.Random.AbstractRNG} <: AbstractMultiplicationMatricesSolver
    ɛ::T
    rng::RNGT
end
ReorderedSchurMultiplicationMatricesSolver(ɛ) = ReorderedSchurMultiplicationMatricesSolver(ɛ, Compat.Random.GLOBAL_RNG)
ReorderedSchurMultiplicationMatricesSolver{T}() where T = ReorderedSchurMultiplicationMatricesSolver(Base.rtoldefault(real(T)))

function solvemultiplicationmatrices(Ms::AbstractVector{<:AbstractMatrix{T}}, solver::ReorderedSchurMultiplicationMatricesSolver) where T
    λ = rand(solver.rng, length(Ms))
    λ /= sum(λ)
    _solvemultiplicationmatrices(Ms, λ, solver)
end

# Deterministic part
function _solvemultiplicationmatrices(Ms::AbstractVector{<:AbstractMatrix{T}}, λ, solver::ReorderedSchurMultiplicationMatricesSolver) where T<:Real
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
    vals
end

function algebraicsolver(algo::AbstractMultiplicationMatricesAlgorithm,
                         solver::AbstractMultiplicationMatricesSolver)
    SolverUsingMultiplicationMatrices(algo, solver)
end

defaultmultiplicationmatricesalgorithm(p) = GröbnerBasisMultiplicationMatricesAlgorithm()
defaultmultiplicationmatricessolver(::Type{T}) where T = ReorderedSchurMultiplicationMatricesSolver{T}()
defaultmultiplicationmatricessolver(::AbstractVector{PT}) where {T, PT<:APL{T}} = defaultmultiplicationmatricessolver(T)

function defaultalgebraicsolver(p)
    algebraicsolver(defaultmultiplicationmatricesalgorithm(p), defaultmultiplicationmatricessolver(p))
end
function defaultalgebraicsolver(p, algo::AbstractMultiplicationMatricesAlgorithm)
    algebraicsolver(algo, defaultmultiplicationmatricessolver(p))
end
function defaultalgebraicsolver(p, solver::AbstractMultiplicationMatricesSolver)
    algebraicsolver(defaultmultiplicationmatricesalgorithm(p), solver)
end
