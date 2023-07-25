import SemialgebraicSets as SS

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
schur_solver = SS.ReorderedSchurMultiplicationMatricesSolver(
    sqrt(eps(Float64)),
    MersenneTwister(0),
)
newton_solver = SS.NewtonTypeDiagonalization()
