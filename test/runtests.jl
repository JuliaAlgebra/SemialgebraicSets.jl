using Test

using MultivariatePolynomials
const MP = MultivariatePolynomials

using SemialgebraicSets

# Taken from JuMP/test/solvers.jl
function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

if try_import(:DynamicPolynomials)
    Mod = DynamicPolynomials
    include("commutativetests.jl")
end

if try_import(:TypedPolynomials)
    Mod = TypedPolynomials
    include("commutativetests.jl")
end
