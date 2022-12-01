using Test

using MultivariatePolynomials
const MP = MultivariatePolynomials

using SemialgebraicSets

struct DummyAlgo <: SemialgebraicSets.Abstractgröbner_basisAlgorithm end
function SemialgebraicSets.gröbner_basis!(::AbstractVector{<:MP.APL}, ::DummyAlgo)
    error("Dummy algo")
end

struct DummySolver <: SemialgebraicSets.AbstractAlgebraicSolver end
function SemialgebraicSets.solve(
    ::SemialgebraicSets.AlgebraicSet,
    ::DummySolver,
)
    error("Dummy solver")
end

function algo_and_solver_test(func)
    Mod.@polyvar x y
    p = [y^2 - x * y, x^2 + x * y]
    V = func(p, DummyAlgo())
    @test_throws ErrorException("Dummy algo") p[1] in V
    @test_throws ErrorException("Dummy algo") rem(p[1], V.I)
    V = func(p, DummySolver())
    @test_throws ErrorException("Dummy solver") collect(V)
    V = func(p, DummyAlgo(), DummySolver())
    @test_throws ErrorException("Dummy solver") collect(V)
    @test_throws ErrorException("Dummy algo") rem(p[1], V.I)
end

@testset "Algo and solver" begin
    algo_and_solver_test(algebraic_set)
    algo_and_solver_test(projective_algebraic_set)
end
