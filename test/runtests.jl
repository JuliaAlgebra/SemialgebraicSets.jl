using SemialgebraicSets
using Base.Test

using DynamicPolynomials
using MultivariatePolynomials

@testset "Basic semialgebraic set" begin
    @polyvar x y
    @test isa(FullSpace(), FullSpace)
    V = @set x * y == 1
    @test V isa AlgebraicSet
    @test_throws ArgumentError addinequality!(V, x*y)
    S = @set x == 0 && x^2*y >= 1
    addequality!(S, x^2 - y)
    addinequality!(S, x + y - 1)
    @test S isa BasicSemialgebraicSet
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

# CLO15
# Ideals, Varieties, and Algorithms
# Cox, Little and O'Shea, Fourth edition, 2015

# CGT97
# Corless, R. M.; Gianni, P. M. & Trager, B. M.
# A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots
# Proceedings of the 1997 international symposium on Symbolic and algebraic computation, 1997, 133-140

# MD95
# Manocha, D. & Demmel, J. Algorithms for intersecting parametric and algebraic curves II: multiple intersections
# Graphical Models and Image Processing, Elsevier, 1995, 57, 81-100


# Taken from CLO15
# They have been adapted to the grlex ordering

#@testset "Equality between algebraic sets" begin
#    @polyvar x y
#    @test (@set x^2 == y && y + x^2 == 4) == (@set x^2 == y && x^2 == 2)
#    @test (@set y == x^2 && x*z == y^2) == (@set y == x^2 && x*z == x^4)
#end

@testset "S-polynomial" begin
    @polyvar x y z
    @test spolynomial(x^3*y^2 - x^2*y^3 + x, 3x^4*y + y^2) == -x^3*y^3 + x^2 - y^3/3
    @test spolynomial(y - x^2, z - x^3) == z - x*y
    @test spolynomial(x^3 - 2x*y, x^2*y - 2y^2 + x) == -x^2
end

@testset "Groebner basis" begin
    @polyvar x y z
    function testg(G, H)
        presort!(G)
        _isz(f) = isapproxzero(rem(f, H))
        @test all(_isz.(G - H))
    end
    function testf(F, H)
        testg(gröbnerbasis(F), H)
    end
    testf([4x^2 + 5x, 3x^3], [x])
    testf([y - x^2, z - x^3], [x*z - y^2, x*y - z, x^2 - y, y^3 - z^2])
    testf([x^3 - 2x*y, x^2*y - 2y^2 + x], [y^2 - 0.5x, x*y, x^2])
    #p = 9x^11 + 27x^14 + x
    # root x=>-0.83005, y=>-3*(-0.83005)^4)
    F = [x^3*y^2 - x^2*y^3 + x, 3x^4*y + y^2]
    H = [x*y^3 - y^4 - 3x^3,
         x^3*y^2 - x^2*y^3 + x,
         x^4*y + y^2/3,
         x^5 + x*y/3,
         y^6 + 3y^5 + 9x^4 + 9x^3*y + y^3/3 - x^2 - x*y - y^2 - 3x]
    function testb(pre, sel)
        testg(groebnerbasis(F, Buchberger(pre, sel)), H)
    end
    testb(identity, dummyselection)
    testb(identity, normalselection)
    testb(presort!, dummyselection)
    testb(presort!, normalselection)
end

function testelements(X, Y; atol=Base.rtoldefault(Float64), kwargs...)
    @test length(X) == length(Y)
    for x in X
        found = false
        for y in Y
            if isapprox(x, y; atol=atol, kwargs...)
                found = true
                break
            end
        end
        @test found
    end
    for y in Y
        found = false
        for x in X
            if isapprox(x, y; atol=atol, kwargs...)
                found = true
                break
            end
        end
        @test found
    end
end

@testset "Example 5.1 of CGT97" begin
    ɛ = 1e-4
    Iɛ = [1 - ɛ 0
          0     1 + ɛ]
    J = [0 1
         1 0]
    Z = zeros(2, 2)
    A = [Iɛ Z
         Z  J]
    B = [J Z
         Z Iɛ]
    α = 0.219
    testelements(SemialgebraicSets._solvemultiplicationmatrices([A, B], [α, 1-α], ReorderedSchurMultiplicationMatricesSolver()), [[1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]; rtol=1e-7)
end

@testset "Example 5.2 of CGT97" begin
    @polyvar x y z
    V = @set x^2 + y^2 == 1 && x^3 + (2 + z)*x*y + y^3 == 1 && z^2 == 2
    @test iszerodimensional(V)
    testelements(V, [[0, 1, √2], [0, 1, -√2], [1, 0, -√2], [1, 0, √2], [-√2/2, -√2/2, √2], [√2/2, √2/2, -√2]])
end

@testset "Example 4.3 of MD95" begin
    @polyvar x y
    V = @set x^2 + 4y^4 == 4y^2 && 3x^4 + y^2 == 5x^3
    # This test is tricky because in the schur decomposition, the 4 last eigenvalues are e.g. 3.4e-7, -1.7e-7+3e-7im, -1.7e-7-3e-7im, -6e-16
    # the second and third do not seem that close but when the three first are averaged it is very close to zero.
    V.solver = algebraicsolver(ReorderedSchurMultiplicationMatricesSolver(1e-5))
    @test iszerodimensional(V)
    testelements(V, [[0.66209555, 0.935259169], [0.66209555, -0.935259169], [0.0516329456, -0.025825086], [0.0516329456, 0.025825086], [0, 0]])
end

@testset "Example 5.3 of CGT97" begin
    @polyvar x y z
    V = @set x^2 + y^2 == 1 && x^3 + (2 + z)*x*y + y^3 == 1 && z^2 == 2
    @test iszerodimensional(V)
    testelements(V, [[0, 1, √2], [0, 1, -√2], [1, 0, -√2], [1, 0, √2], [-√2/2, -√2/2, √2], [√2/2, √2/2, -√2]])
end

@testset "Zero-dimensional ideal" begin
    @polyvar x y z
    V = @set x == y
    @test !iszerodimensional(V)
    @test_throws ErrorException start(V)
    @test_throws ErrorException length(V)
    V =  @set 4x^2 == -5x && 3x^3 == 0
    @test iszerodimensional(V)
    testelements(V, [[0]])
    V = @set y == x^2 && z == x^3
    @test !iszerodimensional(V)
    V = @set x^3 == 2x*y &&  x^2*y == 2y^2 + x
    @test iszerodimensional(V)
    testelements(V, [[0, 0]])
    V = @set x == 1
    @test iszerodimensional(V)
    testelements(V, [[1]])
    V = @set x == 1 && y == 2
    @test iszerodimensional(V)
    testelements(V, [[1, 2]])
    V = @set x == 4 && y^2 == x
    @test iszerodimensional(V)
    testelements(V, [[4, 2], [4, -2]])
    V = @set x^2 + x == 6 && y == x+1
    @test iszerodimensional(V)
    testelements(V, [[2, 3], [-3, -2]])
    V = @set x^2 + x == 6 && y^2 == x
    @test iszerodimensional(V)
    testelements(V, [[2, √2], [2, -√2]])
end
