using SemialgebraicSets
using Base.Test

using DynamicPolynomials

@testset "Basic semialgebraic set" begin
    @polyvar x y
    @test isa(FullSpace(), FullSpace)
    V = @set x * y == 1
    @test V isa AlgebraicSet{Int}
    @test_throws ArgumentError addinequality!(V, x*y)
    S = @set x == 0 && x^2*y >= 1
    addequality!(S, x^2 - y)
    addinequality!(S, x + y - 1)
    @test S isa BasicSemialgebraicSet{Int}
    @test Int32(2)*x^2*y isa MultivariatePolynomials.AbstractTerm{Int32}
    S = (@set Int32(2)*x^2*y == 0 && 1.0x^2*y >= 0 && (6//3)*x^2*y == -y && 1.5x+y >= 0)
    S2 = S ∩ V
    S3 = V ∩ S
    @test S2.p == S3.p == S.p
    @test S2.V.I.p == S3.V.I.p
    T = (@set x*y^2 == -1 && x^2 + y^2 <= 1)
    V2 = @set T.V && V && x + y == 2.0
    @test V2 isa AlgebraicSet
    @test V2.I.p == [x + y - 2.0; T.V.I.p; V.I.p]
    S4 = @set S && T
    @test S4.p == [S.p; T.p]
    @test S4.V.I.p == [S.V.I.p; T.V.I.p]
end

# Taken from
# Ideals, Varieties, and Algorithms
# Cox, Little and O'Shea, Fourth edition
# They have been adapted to the grlex ordering

#@testset "Equality between algebraic sets" begin
#    @test (@set x^2 == y && y + x^2 == 4) == (@set x^2 == y && x^2 == 2)
#    @test (@set y == x^2 && x*z == y^2) == (@set y == x^2 && x*z == x^4)
#end

@testset "S-polynomial" begin
    @polyvar x y z
    @test spolynomial(x^3*y^2 - x^2*y^3 + x, 3x^4*y + y^2) == -x^3*y^3 + x^2 - y^3/3
    @test spolynomial(y - x^2, z - x^3) == z - x*y
    @test spolynomial(x^3 - 2x*y, x^2*y - 2y^2 + x) == -x^2
end

@testset "S-polynomial" begin
    @polyvar x y z
    @test spolynomial(x^3*y^2 - x^2*y^3 + x, 3x^4*y + y^2) == -x^3*y^3 + x^2 - y^3/3
    @test spolynomial(y - x^2, z - x^3) == z - x*y
    @test spolynomial(x^3 - 2x*y, x^2*y - 2y^2 + x) == -x^2
end
