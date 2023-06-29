# CGT97
# Corless, R. M.; Gianni, P. M. & Trager, B. M.
# A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots
# Proceedings of the 1997 international symposium on Symbolic and algebraic computation, 1997, 133-140

# MD95
# Manocha, D. & Demmel, J. Algorithms for intersecting parametric and algebraic curves II: multiple intersections
# Graphical Models and Image Processing, Elsevier, 1995, 57, 81-100

using LinearAlgebra # for I
@testset "Section 4.1 MD95" begin
    η = 1e-10
    A = [zeros(10, 1) Matrix(I, 10, 10); zeros(1, 10) 0.5]
    A[10, 11] = 0
    A[10, 1] = η
    @test sort.(SemialgebraicSets.clusterordschur(A, sqrt(eps(Float64)))[2]) ==
          [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [11]]
end

@testset "Example 4.1 MD95" begin
    A = [
        0 0 0 1 0 0
        0 0 0 0 1 0
        0 0 0 0 0 1
        -1 0 1 -1 0 0
        -1 0 1 -2 -1 0
        -1 0 1 -2 -2 -1
    ]
    clusters = SemialgebraicSets.clusterordschur(A, sqrt(eps(Float64)))[2]
    @test sort(sort.(clusters); by = first) == [[1, 5, 6], [2]]
end

function testelements(X, Y; atol = Base.rtoldefault(Float64), kwargs...)
    @test length(X) == length(Y)
    for y in Y
        @test any(x -> isapprox(x, y; atol = atol, kwargs...), X)
    end
end
function testelementstypes(X, Y; kwargs...)
    testelements(X, Y; kwargs...)
    for T in [Rational{Int}, Float64]
        if X isa FixedVariablesSet
            U = T
        else
            U = float(T)
        end
        V = similar(X, T)
        testelements(V, Y; kwargs...)
        @test eltype(V) == Vector{U}
        @test collect(V) isa Vector{Vector{U}}
    end
end

# We use a fixed RNG in the tests to decrease nondeterminism. There is still nondeterminism in LAPACK though
using Random
schur_solver = ReorderedSchurMultiplicationMatricesSolver(
    sqrt(eps(Float64)),
    MersenneTwister(0),
)
newton_solver = NewtonTypeDiagonalization()

function zero_dimensional_ideal(solver)
    Mod.@polyvar x y z

    V = @set x == y
    @test !is_zero_dimensional(V)
    @test_throws ErrorException iterate(V)
    @test_throws ErrorException length(V)

    V = @set 4x^2 == -5x && 3x^3 == 0 solver
    @test V.solver.solver === solver
    @test is_zero_dimensional(V)
    testelementstypes(V, [[0]])

    V = @set y == x^2 && z == x^3 solver
    @test !is_zero_dimensional(V)

    if !(solver isa NewtonTypeDiagonalization)
        V = @set x^3 == 2x * y && x^2 * y == 2y^2 + x solver
        @test is_zero_dimensional(V)
        testelementstypes(V, [[0, 0]])
    end

    V = @set x == 1 solver
    @test is_zero_dimensional(V)
    testelementstypes(V, [[1]])

    V = @set x == 1 && y == 2 solver
    @test is_zero_dimensional(V)
    testelementstypes(V, [[1, 2]])

    V = @set x == 4 && y^2 == x solver
    @test is_zero_dimensional(V)
    testelementstypes(V, [[4, 2], [4, -2]])

    V = @set x^2 + x == 6 && y == x + 1 solver
    @test is_zero_dimensional(V)
    testelements(V, [[2, 3], [-3, -2]])

    if !(solver isa NewtonTypeDiagonalization)
        V = @set x^2 + x == 6 && y^2 == x solver
        @test is_zero_dimensional(V)
        testelements(V, [[2, √2], [2, -√2]])
    end
end

@testset "Zero-dimensional ideal" begin
    zero_dimensional_ideal(schur_solver)
    zero_dimensional_ideal(newton_solver)
end

function projective_zero_dimensional_ideal(solver)
    Mod.@polyvar x y z

    V = projective_algebraic_set([x - y], solver)
    @test is_zero_dimensional(V)
    testelementstypes(V, [[1, 1]])

    V = @set x + y == z solver
    V.projective = true
    @test !is_zero_dimensional(V)

    V = @set y == 2x solver
    V.projective = true
    @test is_zero_dimensional(V)
    testelementstypes(V, [[1, 2]])

    V = @set x + y == y solver
    V.projective = true
    @test is_zero_dimensional(V)
    testelementstypes(V, [[0, 1]])

    V = projective_algebraic_set([x + y - x])
    @test is_zero_dimensional(V)
    testelementstypes(V, [[1, 0]])
    return
end

@testset "Projective zero-dimensional ideal" begin
    projective_zero_dimensional_ideal(schur_solver)
    projective_zero_dimensional_ideal(newton_solver)
end

function cgt96_e51(solver)
    ɛ = 1e-4
    Iɛ = [
        1-ɛ 0
        0 1+ɛ
    ]
    J = [
        0 1
        1 0
    ]
    Z = zeros(2, 2)
    A = [
        Iɛ Z
        Z J
    ]
    B = [
        J Z
        Z Iɛ
    ]
    α = 0.219
    return testelements(
        SemialgebraicSets._solve_multiplication_matrices(
            [A, B],
            [α, 1 - α],
            solver,
        ),
        [[1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
        rtol = 1e-7,
    )
end

@testset "Example 5.1 of CGT97" begin
    @testset "Schur" begin
        cgt96_e51(schur_solver)
    end
    #    @testset "Newton" begin
    #        cgt96_e51(newton_solver)
    #    end
end

function md95_e43(solver)
    Mod.@polyvar x y
    V = @set x^2 + 4y^4 == 4y^2 && 3x^4 + y^2 == 5x^3 solver
    # This test is tricky because in the schur decomposition, the 4 last eigenvalues are e.g. 3.4e-7, -1.7e-7+3e-7im, -1.7e-7-3e-7im, -6e-16
    # the second and third do not seem that close but when the three first are averaged it is very close to zero.
    @test is_zero_dimensional(V)
    return testelementstypes(
        V,
        [
            [0.66209555, 0.935259169],
            [0.66209555, -0.935259169],
            [0.0516329456, -0.025825086],
            [0.0516329456, 0.025825086],
            [0, 0],
        ],
    )
end

@testset "Example 4.3 of MD95" begin
    md95_e43(schur_solver)
end

@testset "Example 5.2 of CGT97" begin
    Mod.@polyvar x y z
    V =
        @set x^2 + y^2 == 1 && x^3 + (2 + z) * x * y + y^3 == 1 && z^2 == 2 schur_solver
    @test is_zero_dimensional(V)
    B = standard_monomials(V.I)
    @test !isnothing(B)
    @test B == [
        y^3 * z,
        x * y * z,
        y^3,
        y^2 * z,
        x * y,
        x * z,
        y^2,
        y * z,
        x,
        y,
        z,
        1,
    ]
    testelements(
        V,
        [
            [0, 1, √2],
            [0, 1, -√2],
            [1, 0, -√2],
            [1, 0, √2],
            [-√2 / 2, -√2 / 2, √2],
            [√2 / 2, √2 / 2, -√2],
        ],
    )
end

#@testset "Example 4.4 of MD95 and 5.3 of CGT97" begin
#    Mod.@polyvar x y
#    F = -2-7x+14x^3-7x^5+x^7 + (7-42x^2+35x^4-7x^6)*y + (16+42x-70x^3+21x^5)*y^2 + (-14+70x^2-35x^4)*y^3 + (-20-35x+35x^3)*y^4 + (7-21x^2)*y^5 + (8+7x)*y^6 - y^7 - y^8
#    G = differentiate(F, y)
#    @show G
#    V = @set F == 0 && G == 0 solver
#    @test is_zero_dimensional(V)
#    sols = [[-3.91298, -1.95065], [-3.23984, -1.56367], [-3.21615, -1.41421], [-3.09474, -1.84776], [3.6497, 1.84776], [3.65578, 1.80399], [2.2928, 1.84776], [2.66119, 1.41421], [2.68379, 1.23369], [2.5673, 0.765367], [0.600779, 1.84776], [0.969172, 1.41421], [-1.85925, -1.41422], [-1.40272, -1.84775], [-2.01312, -0.812102], [-2.01235, -0.765367], [1.26106, 0.265359], [-0.481613, 0.765367], [1.03657, -0.765367], [-0.483778, 0.630692], [-0.0458213, -1.84776], [-1.80194, 0.]]
#    @show length(sols)
#    #testelements(V, sols)
#    @show length(V)
#    for v in collect(V)
#        @show v
#        @show F(variables(F) => v)
#        @show G(variables(G) => v)
#    end
#end

@testset "Complex not yet implemented" begin
    Mod.@polyvar x
    V = @set x^2 == im
    @test_throws MethodError collect(V)
end
