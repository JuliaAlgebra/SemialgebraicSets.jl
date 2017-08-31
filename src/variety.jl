export iszerodimensional
export algebraicset, equalities

type DefaultAlgebraicSetLibrary{S<:AbstractAlgebraicSolver}
    solver::S
end

defaultalgebraicsetlibrary(::Vector{<:APL}, solver::AbstractAlgebraicSolver) = DefaultAlgebraicSetLibrary(solver)
defaultalgebraicsetlibrary(p::Vector{<:APL}, solveroralgo...) = defaultalgebraicsolver(p, solveroralgo...)

mutable struct AlgebraicSet{T, PT<:APL{T}, S<:AbstractAlgebraicSolver} <: AbstractAlgebraicSet
    I::PolynomialIdeal{T, PT}
    elements::Vector{Vector{T}}
    elementscomputed::Bool
    iszerodimensional::Bool
    solver::S
end
AlgebraicSet{T, PT, S}(I::PolynomialIdeal{T, PT}, solver::S) where {T, PT, S} = AlgebraicSet{T, PT, S}(I, Vector{T}[], false, false, solver)
AlgebraicSet(I::PolynomialIdeal{T, PT}, solver::S) where {T, PT, S} = AlgebraicSet{T, PT, S}(I, solver)
AlgebraicSet{T, PT}() where {T, PT} = AlgebraicSet(PolynomialIdeal{T, PT}(), defaultalgebraicsolver(T))
AlgebraicSet(p::Vector, solver) = AlgebraicSet(PolynomialIdeal(p), solver)

algebraicset(p::Vector, lib::DefaultAlgebraicSetLibrary) = AlgebraicSet(p, lib.solver)
algebraicset(p::Vector, solver) = algebraicset(p, defaultalgebraicsetlibrary(p, solver))
algebraicset(p::Vector) = algebraicset(p, defaultalgebraicsetlibrary(p))

equalities(V::AlgebraicSet) = V.I.p
addequality!(V::AlgebraicSet, p) = push!(V.I.p, p)
Base.intersect(S::AlgebraicSet, T::AlgebraicSet) = AlgebraicSet(S.I + T.I, S.solver)

defaultalgebraicsolver(V::AlgebraicSet) = V.solver
function computeelements!(V::AlgebraicSet{T}) where T
    if !V.elementscomputed
        elements = solvealgebraicequations(V, V.solver)
        V.iszerodimensional = !isnull(elements)
        if !isnull(elements)
            V.elements = get(elements)
        end
        V.elementscomputed = true
    end
end
function iszerodimensional(V::AlgebraicSet)
    computeelements!(V)
    V.iszerodimensional
end

Base.eltype(V::AlgebraicSet{T}) where T = Vector{T}

for f in [:length, :start, :endof, :next, :done, :getindex]
    @eval begin
        function Base.$f(V::AlgebraicSet, args...)
            computeelements!(V)
            if !iszerodimensional(V)
                error("A non zero-dimensional algebraic set is not iterable")
            end
            $f(V.elements, args...)
        end
    end
end
