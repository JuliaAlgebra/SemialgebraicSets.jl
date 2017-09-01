@testset "Basic semialgebraic set" begin
    Mod.@polyvar x y
    @test isa(FullSpace(), FullSpace)
    V = @set x * y == 1
    @test V isa AlgebraicSet
    @test_throws ArgumentError addinequality!(V, x*y)
    S = @set x == 0 && x^2*y >= 1
    addequality!(S, x^2 - y)
    addinequality!(S, x + y - 1)
    @test S isa BasicSemialgebraicSet
    @test BasicSemialgebraicSet{Int, polynomialtype(x, Int)}() isa BasicSemialgebraicSet{Float64, polynomialtype(x, Float64)}
    @test Int32(2)*x^2*y isa MultivariatePolynomials.AbstractTerm{Int32}
    S = (@set Int32(2)*x^2*y == 0 && 1.0x^2*y >= 0 && (6//3)*x^2*y == -y && 1.5x+y >= 0)
    S2 = S ∩ V
    S3 = V ∩ S
    @test inequalities(S2) == inequalities(S3) == S.p
    @test equalities(S2) == S3.V.I.p
    T = (@set x*y^2 == -1 && x^2 + y^2 <= 1)
    V2 = @set T.V && V && x + y == 2.0
    @test V2 isa AlgebraicSet
    @test V2.I.p == [x + y - 2.0; equalities(T); equalities(V)]
    S4 = @set S && T
    @test S4.p == [S.p; inequalities(T)]
    @test equalities(S4) == [S.V.I.p; T.V.I.p]
end
