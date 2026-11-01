export ReorderedSchurMultiplicationMatricesSolver

struct MultiplicationMatrices{Ms}
    matrices::Ms
end

"""
    AbstractMultiplicationMatricesSolver

Solver of algebraic equations using multiplication matrices.
"""
abstract type AbstractMultiplicationMatricesSolver end

function solve(
    Ms::MultiplicationMatrices,
    solver::AbstractMultiplicationMatricesSolver,
)
    λ = rand(solver.rng, length(Ms.matrices))
    λ /= sum(λ)
    return _solve_multiplication_matrices(Ms.matrices, λ, solver)
end

include("cluster.jl")
include("schur.jl")
include("newton_type.jl")

function default_multiplication_matrices_solver(
    ::AbstractVector{PT},
) where {PT<:_APL}
    return default_multiplication_matrices_solver(MP.coefficient_type(PT))
end

function default_multiplication_matrices_solver(::Type{T}) where {T}
    return ReorderedSchurMultiplicationMatricesSolver{T}()
end
