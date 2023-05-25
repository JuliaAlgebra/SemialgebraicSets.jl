using Test

import MultivariatePolynomials as MP

using SemialgebraicSets

struct DummyAlgo <: SemialgebraicSets.AbstractGröbnerBasisAlgorithm end
SemialgebraicSets.default_algebraic_solver(::Any, ::DummyAlgo) = DummySolver()
SemialgebraicSets.promote_for(::Type, ::Type{DummyAlgo}) = Int
function SemialgebraicSets.gröbner_basis!(
    ::AbstractVector{<:MP.APL},
    ::DummyAlgo,
)
    return error("Dummy algo")
end

struct DummySolver <: SemialgebraicSets.AbstractAlgebraicSolver end
function SemialgebraicSets.default_gröbner_basis_algorithm(::Any, ::DummySolver)
    return SemialgebraicSets.NoAlgorithm()
end
SemialgebraicSets.promote_for(::Type{T}, ::Type{DummySolver}) where {T} = T
function SemialgebraicSets.solve(
    ::SemialgebraicSets.AlgebraicSet,
    ::DummySolver,
)
    return error("Dummy solver")
end

function algo_and_solver_test(func)
    Mod.@polyvar x y
    p = [y^2 - x * y, x^2 + x * y]
    V = func(p, DummyAlgo(), SemialgebraicSets.default_algebraic_set_library(p))
    @test V.I.algo isa DummyAlgo
    @test_throws ErrorException("Dummy algo") p[1] in V
    @test_throws ErrorException("Dummy algo") rem(p[1], V.I)
    V = func(p, DummyAlgo())
    @test V.I.algo isa DummyAlgo
    @test V.solver isa DummySolver
    @test eltype(V) == Vector{Int}
    @test_throws ErrorException("Dummy solver") p[1] in V
    @test_throws ErrorException("Dummy algo") rem(p[1], V.I)
    V = func(p, DummySolver())
    @test V.I.algo isa SemialgebraicSets.NoAlgorithm
    @test V.solver isa DummySolver
    @test eltype(V) == Vector{Int}
    @test_throws ErrorException("Dummy solver") collect(V)
    V = func(p, DummyAlgo(), DummySolver())
    @test V.I.algo isa DummyAlgo
    @test V.solver isa DummySolver
    @test eltype(V) == Vector{Int}
    @test_throws ErrorException("Dummy solver") collect(V)
    @test_throws ErrorException("Dummy algo") rem(p[1], V.I)
end

@testset "Algo and solver" begin
    algo_and_solver_test(algebraic_set)
    #algo_and_solver_test(projective_algebraic_set)
end
