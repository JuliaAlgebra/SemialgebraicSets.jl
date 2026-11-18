using Test
using SemialgebraicSets

import MultivariatePolynomials as MP
import DynamicPolynomials

Mod = isdefined(Main, :Mod) ? Main.Mod : DynamicPolynomials

@testset "Truncated" begin
    Mod.@polyvar x y
    V = @set x * y == 1 && y^2 == 1
    W = truncated(V, 3)
    @test W isa SemialgebraicSets.TruncatedAlgebraicSet
    @test W isa AbstractAlgebraicSet
    @test equalities(W) == equalities(V)
    @test SemialgebraicSets.nequalities(W) == 2
    @test MP.variables(W) == MP.variables(V)
    @test MP.monomial_type(typeof(W)) == MP.monomial_type(typeof(V))
    @test MP.coefficient_type(typeof(W)) == MP.coefficient_type(typeof(V))
    I = ideal(W)
    @test I isa SemialgebraicSets.TruncatedIdeal
    @test I.generators == equalities(V)
    @test I.maxdegree == 3
    @test MP.coefficient_type(typeof(I)) == MP.coefficient_type(typeof(V))
    # The Gröbner basis is not computed
    @test !SemialgebraicSets.is_gröbner_basis(ideal(V))
    W2 = similar(W, Float64)
    @test W2 isa SemialgebraicSets.TruncatedAlgebraicSet{Float64}
    @test W2.maxdegree == 3
    @test similar(W2, Float64) === W2
    @test sprint(show, W) == "truncated($(sprint(show, V)), 3)"
    @test occursin("truncated at degree 3", sprint(show, MIME"text/plain"(), W))
end
