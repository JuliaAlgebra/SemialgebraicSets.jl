using Test

import StarAlgebras as SA
import MultivariatePolynomials as MP
import MultivariateBases as MB
import DynamicPolynomials

using SemialgebraicSets

Mod = isdefined(Main, :Mod) ? Main.Mod : DynamicPolynomials

@testset "promote_bases_with_maps" begin
    Mod.@polyvar x y z

    @testset "FullSpace with polynomial" begin
        f = FullSpace()
        p = x^2 + y
        (new_f, map_f), (new_p, map_p) = SA.promote_bases_with_maps(f, p)
        @test new_f isa FullSpace
        @test map_f === nothing
        @test new_p === p
        @test map_p === nothing
    end

    @testset "FullSpace with AlgebraicSet" begin
        f = FullSpace()
        V = @set x * y == 1
        (new_f, map_f), (new_V, map_V) = SA.promote_bases_with_maps(f, V)
        @test new_f isa FullSpace
        @test map_f === nothing
        @test new_V === V
        @test map_V === nothing
    end

    @testset "FullSpace with FullSpace" begin
        f1 = FullSpace()
        f2 = FullSpace()
        (new_f1, map_f1), (new_f2, map_f2) = SA.promote_bases_with_maps(f1, f2)
        @test new_f1 isa FullSpace
        @test map_f1 === nothing
        @test new_f2 isa FullSpace
        @test map_f2 === nothing
    end

    @testset "AlgebraicSet with polynomial - same variables" begin
        V = @set x * y == 1
        p = x + y
        (new_V, map_V), (new_p, map_p) = SA.promote_bases_with_maps(V, p)
        @test new_V === V
        @test map_V === nothing
        @test new_p === p
        @test map_p === nothing
    end

    @testset "AlgebraicSet with polynomial - different variables" begin
        V = @set x^2 == 1
        p = y + z
        (new_V, map_V), (new_p, map_p) = SA.promote_bases_with_maps(V, p)
        @test map_V !== nothing
        @test map_p !== nothing
        @test Set(MP.variables(new_V)) == Set([x, y, z])
        new_eq = equalities(new_V)
        @test length(new_eq) == 1
        @test Set(MP.variables(new_eq[1])) == Set([x, y, z])
        @test new_eq[1] == x^2 - 1
        @test new_p == p
        @test Set(MP.variables(p)) == Set([y, z])
    end

    @testset "AlgebraicSet with AlgebraicSet - same variables" begin
        V1 = @set x * y == 1
        V2 = @set x + y == 2
        (new_V1, map_V1), (new_V2, map_V2) = SA.promote_bases_with_maps(V1, V2)
        @test new_V1 === V1
        @test map_V1 === nothing
        @test new_V2 === V2
        @test map_V2 === nothing
    end

    @testset "AlgebraicSet with term" begin
        V = @set x^2 == 1
        p = 2y
        (new_V, map_V), (new_p, map_p) = SA.promote_bases_with_maps(V, p)
        @test Set(MP.variables(new_p)) == Set([x, y])
        @test new_p == p
        @test only(equalities(new_V)) == x^2 - 1
        @test map_V !== nothing
        @test map_p !== nothing
    end

    @testset "AlgebraicSet with AlgebraicSet - different variables" begin
        V1 = @set x^2 == 1
        V2 = @set y^2 == 4
        (new_V1, map_V1), (new_V2, map_V2) = SA.promote_bases_with_maps(V1, V2)
        @test map_V1 !== nothing
        @test map_V2 !== nothing
        @test Set(MP.variables(new_V1)) == Set([x, y])
        @test Set(MP.variables(new_V2)) == Set([x, y])
        eq1 = equalities(new_V1)
        @test length(eq1) == 1
        @test Set(MP.variables(eq1[1])) == Set([x, y])
        eq2 = equalities(new_V2)
        @test length(eq2) == 1
        @test Set(MP.variables(eq2[1])) == Set([x, y])
    end

    @testset "BasicSemialgebraicSet with polynomial - same variables" begin
        S = @set x - y == 0 && x^2 * y >= 1
        p = x + y
        (new_S, map_S), (new_p, map_p) = SA.promote_bases_with_maps(S, p)
        @test new_S === S
        @test map_S === nothing
        @test new_p === p
        @test map_p === nothing
    end

    @testset "BasicSemialgebraicSet with polynomial - different variables" begin
        S = @set x >= 1
        p = y + z
        (new_S, map_S), (new_p, map_p) = SA.promote_bases_with_maps(S, p)
        @test map_S !== nothing
        @test map_p !== nothing
        ineqs = inequalities(new_S)
        @test length(ineqs) == 1
        @test Set(MP.variables(ineqs[1])) == Set([x, y, z])
    end

    @testset "BasicSemialgebraicSet with AlgebraicSet - same variables" begin
        S = @set x - y == 0 && x^2 * y >= 1
        V = @set x + y == 2
        (new_S, map_S), (new_V, map_V) = SA.promote_bases_with_maps(S, V)
        @test new_S === S
        @test map_S === nothing
        @test new_V === V
        @test map_V === nothing
    end

    @testset "BasicSemialgebraicSet with BasicSemialgebraicSet - same variables" begin
        S1 = @set x - y == 0 && x^2 * y >= 1
        S2 = @set x + y == 2 && x * y >= 0
        (new_S1, map_S1), (new_S2, map_S2) = SA.promote_bases_with_maps(S1, S2)
        @test new_S1 === S1
        @test map_S1 === nothing
        @test new_S2 === S2
        @test map_S2 === nothing
    end

    @testset "BasicSemialgebraicSet with BasicSemialgebraicSet - different variables" begin
        S1 = @set x >= 1
        S2 = @set y >= 2
        (new_S1, map_S1), (new_S2, map_S2) = SA.promote_bases_with_maps(S1, S2)
        @test map_S1 !== nothing
        @test map_S2 !== nothing
        ineqs1 = inequalities(new_S1)
        @test length(ineqs1) == 1
        @test Set(MP.variables(ineqs1[1])) == Set([x, y])
        ineqs2 = inequalities(new_S2)
        @test length(ineqs2) == 1
        @test Set(MP.variables(ineqs2[1])) == Set([x, y])
    end

    @testset "BasicSemialgebraicSet (with FullSpace V) - different variables" begin
        S = @set x >= 1
        p = y^2
        @test S.V isa FullSpace
        (new_S, map_S), (new_p, map_p) = SA.promote_bases_with_maps(S, p)
        @test map_S !== nothing
        @test map_p !== nothing
        @test new_S.V isa FullSpace
        ineqs = inequalities(new_S)
        @test length(ineqs) == 1
        @test Set(MP.variables(ineqs[1])) == Set([x, y])
    end

    @testset "FixedVariablesSet with polynomial - different variables" begin
        V = FixedVariablesSet(
            SemialgebraicSets.FixedVariablesIdeal{typeof(x),Int,typeof(x * y)}(
                Dict(x => 1),
            ),
        )
        p = y + z
        (new_V, map_V), (new_p, map_p) = SA.promote_bases_with_maps(V, p)
        @test map_V !== nothing
        @test map_p !== nothing
        @test new_V === V
    end

    @testset "Polynomials promoted to common variables at construction" begin
        S = @set x >= 1 && y >= 1
        ineqs = inequalities(S)
        @test MP.variables(ineqs[1]) == MP.variables(ineqs[2])
        @test Set(MP.variables(ineqs[1])) == Set([x, y])
        @test ineqs == [x - 1, y - 1]

        V = @set x^2 == 1 && y^2 == 4
        eqs = equalities(V)
        @test MP.variables(eqs[1]) == MP.variables(eqs[2])
        @test Set(MP.variables(eqs[1])) == Set([x, y])
    end

    @testset "AlgebraicSet with polynomial - overlapping variables" begin
        V = @set x + y == 1
        p = y + z
        (new_V, map_V), (new_p, map_p) = SA.promote_bases_with_maps(V, p)
        @test map_V !== nothing
        @test map_p !== nothing
        eq = equalities(new_V)
        @test length(eq) == 1
        @test Set(MP.variables(eq[1])) == Set([x, y, z])
    end

    @testset "Preserve $B basis" for B in (
        MP.Monomial,
        MB.ScaledMonomial,
        MB.Chebyshev,
    )
        full = MP.FullBasis{B}(MP.variables(y * z))
        index = MP.exponents(y * z^2)
        keys = [zero.(index), index]
        for basis in (full, SA.SubBasis(full, keys))
            coefs =
                basis === full ? SA.SparseCoefficients(keys, [1.0, 2.0]) :
                [1.0, 2.0]
            p = SA.AlgebraElement(coefs, MP.algebra(basis))
            index = basis === full ? keys[2] : 2
            t = SA.Term(parent(p), index, 2.0)
            for operand in (p, t)
                (space, map_space), (same, map_same) =
                    SA.promote_bases_with_maps(FullSpace(), operand)
                @test space isa FullSpace
                @test map_space === nothing
                @test same === operand
                @test map_same === nothing

                domain = @set y + z >= 0
                (same_domain, map_domain), (same, map_same) =
                    SA.promote_bases_with_maps(domain, operand)
                @test same_domain === domain
                @test map_domain === nothing
                @test same === operand
                @test map_same === nothing

                for domain in
                    (@set(x^2 == 1), @set(x >= 0), @set(x^2 == 1 && x >= 0))
                    (new_domain, map_domain), (q, map_q) =
                        SA.promote_bases_with_maps(domain, operand)
                    @test MP.variables(new_domain) == MP.variables(x * y * z)
                    @test MP.variables(q) == MP.variables(x * y * z)
                    @test MP.variables(operand) == MP.variables(y * z)
                    @test MP.variables(domain) == MP.variables(x)
                    @test q == operand
                    @test eltype(SA.basis(parent(q))) <: MP.BasisElement{B}
                    @test (SA.basis(parent(q)) isa SA.SubBasis) ==
                          (basis isa SA.SubBasis)
                    @test map_domain !== nothing
                    @test map_q !== nothing
                    @test (q isa SA.Term) == (operand isa SA.Term)
                    if operand isa SA.Term
                        @test map_q(operand.index) == q.index
                    else
                        @test map(
                            map_q,
                            collect(SA.keys(SA.coeffs(operand))),
                        ) == collect(SA.keys(SA.coeffs(q)))
                    end
                end
            end
        end
    end

    @testset "Noncommutative domain promotion" begin
        DynamicPolynomials.@ncpolyvar a b c
        domain = @set a >= 0
        p = b * c + c * b
        new_domain, q = SA.promote_bases(domain, p)
        @test MP.variables(new_domain) == [a, b, c]
        @test MP.variables(q) == [a, b, c]
        @test q == p
        @test q isa typeof(p)
        @test MP.variables(p) == [b, c]
    end
end
