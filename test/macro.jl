struct DummySolver <: SemialgebraicSets.AbstractAlgebraicSolver end

@testset "Basic semialgebraic set" begin
    Mod.@polyvar x y
    @test isa(FullSpace(), FullSpace)
    V = @set x * y == 1
    @test V isa AlgebraicSet
    @test_throws ArgumentError addinequality!(V, x*y)
    @testset "Basic" begin
        S = @set x - y == 0 && x^2*y >= 1
        addequality!(S, x^2 - y)
        addinequality!(S, x + y - 1)
        @test S isa BasicSemialgebraicSet
        @test Int32(2)*x^2*y isa MultivariatePolynomials.AbstractTerm{Int32}

        @testset "Mixed types" begin
            S = (@set Int32(2)*x^2*y == 0 && 1.0x^2*y >= 0 && (6//3)*x^2*y == -y && 1.5x+y >= 0)
            S2 = S ∩ V
            S3 = V ∩ S
            @test inequalities(S2) == inequalities(S3) == S.p
            @test equalities(S2) == S3.V.I.p
        end

        T = (@set x*y^2 == -1 && x^2 + y^2 <= 1)
        V2 = @set T.V && V && x + y == 2.0
        @test V2 isa AlgebraicSet
        @test V2.I.p == [equalities(T); equalities(V); x + y - 2.0]
        S4 = @set S && T
        @test S4.p == [S.p; inequalities(T)]
        @test equalities(S4) == [S.V.I.p; T.V.I.p]
    end
    @testset "Basic with no equality" begin
        S = @set x + y ≥ 1 && x ≤ y
        @test S isa BasicSemialgebraicSet
        @test S.V isa FullSpace
    end
    @testset "Basic with fixed variables" begin
        S = @set x == 1 && x ≤ y && 2 == y
        @test S isa BasicSemialgebraicSet
        @test S.V isa FixedVariablesSet
        @test rem(x + y, ideal(S.V)) == 3

        S = @set x == 1 && x ≤ y && 2 == y && 1 == x
        @test S isa BasicSemialgebraicSet
        @test S.V isa FixedVariablesSet
        @test rem(x + y, ideal(S.V)) == 3

        @testset "Empty" begin
            for S in (@set(x == 1 && x ≤ y && 2 == y && 2 == x),
                      @set(x == 1 && x ≤ y && 2 == x),
                      @set(x == 1 && x ≤ y && 2 == x && y == 3),
                      @set(x ≤ y && y == 3 && x == 1 && y == 1),
                      @set(x ≤ y && y == 3 && x == 1 && y == 3 && x == 2))
                @test S isa BasicSemialgebraicSet
                @test S.V isa FixedVariablesSet
                @test isempty(S.V)
                @test iszero(length(S.V))
                @test isempty(collect(S.V))
                @test iszero(rem(x + y, ideal(S.V)))
            end
        end
    end
    @testset "Basic mixed fixed variables and equalities" begin
        for S in [@set(x == 1 && x ≤ y && 2 + y == x),
                  @set(x == 1 && 2 + y == x && x ≤ y),
                  @set(x ≤ y && 2 + y == x && x == 1),
                  @set(x ≤ y && x == 1 && 2 + y == x),
                  @set(2 + y == x && x ≤ y && x == 1),
                  @set(2 + y == x && x ≤ y && x == 1)]
            @test S isa BasicSemialgebraicSet
            @test S.V isa AlgebraicSet
        end

        solver = DummySolver()
        for S in [@set(x == 1 && x ≤ y && 2 + y == x, solver),
                  @set(x == 1 && 2 + y == x && x ≤ y, solver),
                  @set(x ≤ y && 2 + y == x && x == 1, solver),
                  @set(x ≤ y && x == 1 && 2 + y == x, solver),
                  @set(2 + y == x && x ≤ y && x == 1, solver),
                  @set(2 + y == x && x ≤ y && x == 1, solver)]
            @test S isa BasicSemialgebraicSet
            @test S.V isa AlgebraicSet
            @test S.V.solver === solver
        end
    end
end
