# CLO15
# Ideals, Varieties, and Algorithms
# Cox, Little and O'Shea, Fourth edition, 2015

# Taken from CLO15
# They have been adapted to the grlex ordering

#@testset "Equality between algebraic sets" begin
#    @polyvar x y
#    @test (@set x^2 == y && y + x^2 == 4) == (@set x^2 == y && x^2 == 2)
#    @test (@set y == x^2 && x*z == y^2) == (@set y == x^2 && x*z == x^4)
#end

@testset "S-polynomial" begin
    Mod.@polyvar x y z
    @test spolynomial(x^3*y^2 - x^2*y^3 + x, 3x^4*y + y^2) == -x^3*y^3 + x^2 - y^3/3
    @test spolynomial(y - x^2, z - x^3) == z - x*y
    @test spolynomial(x^3 - 2x*y, x^2*y - 2y^2 + x) == -x^2
end

@testset "Groebner basis" begin
    Mod.@polyvar x y z
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
        testg(groebnerbasis(F, Buchberger(Base.rtoldefault(Float64), pre, sel)), H)
    end
    testb(identity, dummyselection)
    testb(identity, normalselection)
    testb(presort!, dummyselection)
    testb(presort!, normalselection)
end

@testset "Reduce" begin
    Mod.@polyvar x y
    V1 = @set x == 1 && y == x^2
    @test rem(x^2 + 3x*y + 2y, ideal(V1)) == 6
    p = x^2 + x
    V2 = FullSpace()
    @test rem(p, ideal(V2)) === p
    # Needs MP v0.1.1
    #V3 = @set x == y^2
    #@test rem(x^2 + 3x*y + 2y + y^4, ideal(V3)) == rem(2y^4 + 3y^3 + 2y, ideal(V3))
end
