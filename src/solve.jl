export algebraic_solver

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
