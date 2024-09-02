# SemialgebraicSets

| **Build Status** |
|:----------------:|
| [![Build Status][build-img]][build-url] |
| [![Codecov branch][codecov-img]][codecov-url] |

Extension of MultivariatePolynomials to semialgebraic sets, i.e. sets defined by inequalities and equalities between polynomials.
The following example shows how to build an [algebraic set/algebraic variety](https://en.wikipedia.org/wiki/Algebraic_variety)
```julia
julia> using SemialgebraicSets, TypedPolynomials

julia> @polyvar x y z;

julia> @set x*y^2 == x*z - y && x*y == z^2 && x == y*z^4
Algebraic Set defined by 3 equalities
 y - x*z + x*y^2 = 0
 -z^2 + x*y = 0
 x - y*z^4 = 0

julia> algebraic_set([x*y^2 - x*z - y, x*y - z^2, x - y*z^4])
Algebraic Set defined by 3 equalities
 -y - x*z + x*y^2 = 0
 -z^2 + x*y = 0
 x - y*z^4 = 0
```

The following example shows how to build an [basic semialgebraic set](http://www.mit.edu/~parrilo/cdc03_workshop/10_positivstellensatz_2003_12_07_02_screen.pdf).
```julia
julia> using SemialgebraicSets, TypedPolynomials

julia> @polyvar x y;

julia> @set x^2 + y^2 <= 1 # Euclidean ball

julia> @set y^2 == x^3 - x && x <= 0 # Cutting the algebraic variety https://en.wikipedia.org/wiki/Algebraic_variety#/media/File:Elliptic_curve2.png
Basic semialgebraic Set defined by 1 equalitty
 x + y^2 - x^3 = 0
1 inequalitty
 -x ≥ 0


julia> basic_semialgebraic_set(algebraic_set([y^2- x^3 - x]), [-x])
Basic semialgebraic Set defined by 1 equalitty
 -x + y^2 - x^3 = 0
1 inequalitty
 -x ≥ 0
```

## Solving systems of algebraic equations

Once the algebraic set has been created, you can check whether it is zero-dimensional and if it is the case, you can get the finite number of elements of the set simply by iterating over it, or by using `collect` to transform it into an array containing the solutions.
```julia
V = @set y == x^2 && z == x^3
is_zero_dimensional(V) # should return false
V = @set x^2 + x == 6 && y == x+1
is_zero_dimensional(V) # should return true
collect(V) # should return [[2, 3], [-3, -2]]
```
The code sample above solves the system of algbraic equations by first
computing a *Gröbner basis* for the system, then the multiplication matrices
and then a Schur decomposition of a random combination of these matrices.
Additionally, SemialgebraicSets defines an interface that can be implemented by
other solvers for these systems as shown in the following subsections.

### Solve with [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/)

You can solve the system with
[homotopy continuation](https://www.juliahomotopycontinuation.org/)
as follows:
```julia
julia> using HomotopyContinuation

julia> solver = SemialgebraicSetsHCSolver(; compile = false)
SemialgebraicSetsHCSolver(; compile = false)

julia> @polyvar x y
(x, y)

julia> V = @set x^2 + x == 6 && y == x+1 solver
Algebraic Set defined by 2 equalities
 x^2 + x - 6.0 = 0
 -x + y - 1.0 = 0

julia> collect(V)
2-element Vector{Vector{Float64}}:
 [2.0, 3.0]
 [-3.0, -2.0]
```

### Solve with [MacaulayLab](http://www.macaulaylab.net/)

You can solve the system with
[MacaulayLab](http://www.macaulaylab.net/) as follows.
First install [MacaulayLab.jl](https://github.com/blegat/MacaulayLab.jl)
and then run the following:
```julia
julia> using DynamicPolynomial, MacaulayLab, SemialgebraicSets

julia> solver = MacaulayLab.Solver()
MacaulayLab.Solver()

julia> V = @set x^2 + x == 6 && y == x + 1 solver
Algebraic Set defined by 2 equalities
 x^2 + x - 6.0 = 0
 -x + y - 1.0 = 0

julia> collect(V)
2-element Vector{Vector{Float64}}:
 [2.0000000000000004, 2.999999999999999]
 [-3.0000000000000004, -2.0000000000000004]
```

[build-img]: https://github.com/JuliaAlgebra/SemialgebraicSets.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/JuliaAlgebra/SemialgebraicSets.jl/actions?query=workflow%3ACI
[codecov-img]: http://codecov.io/github/JuliaAlgebra/SemialgebraicSets.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaAlgebra/SemialgebraicSets.jl?branch=master
