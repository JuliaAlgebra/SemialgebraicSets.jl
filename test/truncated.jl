using Test
using SemialgebraicSets

import MultivariatePolynomials as MP
import StarAlgebras as SA
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

@testset "Truncated with inequalities" begin
    Mod.@polyvar x y z
    V = @set x * y == 1 && y^2 == 1
    K = truncated(V, 3) ∩ @set(x >= 0)
    @test K isa BasicSemialgebraicSet
    @test algebraic_set(K) isa SemialgebraicSets.TruncatedAlgebraicSet
    @test algebraic_set(K).maxdegree == 3
    @test inequalities(K) == [x]
    # Promotion keeps the truncation
    V = @set x == y
    for p in [x + y, z]
        T = truncated(V, 2)
        (new_T, _), _ = SA.promote_bases_with_maps(T, p)
        @test new_T isa SemialgebraicSets.TruncatedAlgebraicSet
        @test new_T.maxdegree == 2
        @test MP.variables(only(equalities(new_T))) == MP.variables(new_T)
        K = T ∩ @set(y >= 1)
        (new_K, _), _ = SA.promote_bases_with_maps(K, p)
        @test algebraic_set(new_K) isa SemialgebraicSets.TruncatedAlgebraicSet
        @test algebraic_set(new_K).maxdegree == 2
        @test MP.variables(only(equalities(new_K))) == MP.variables(new_K)
        @test MP.variables(only(inequalities(new_K))) == MP.variables(new_K)
    end
end
