export iszerodimensional
export algebraicset, projectivealgebraicset, equalities

type DefaultAlgebraicSetLibrary{S<:AbstractAlgebraicSolver}
    solver::S
end

defaultalgebraicsetlibrary(::Vector{<:APL}, solver::AbstractAlgebraicSolver) = DefaultAlgebraicSetLibrary(solver)
defaultalgebraicsetlibrary(p::Vector{<:APL}, solveroralgo...) = defaultalgebraicsolver(p, solveroralgo...)

mutable struct AlgebraicSet{T, PT<:APL{T}, S<:AbstractAlgebraicSolver} <: AbstractAlgebraicSet
    I::PolynomialIdeal{T, PT}
    projective::Bool
    elements::Vector{Vector{T}}
    elementscomputed::Bool
    iszerodimensional::Bool
    solver::S
end
AlgebraicSet{T, PT, S}(I::PolynomialIdeal{T, PT}, solver::S) where {T, PT, S} = AlgebraicSet{T, PT, S}(I, false, Vector{T}[], false, false, solver)
AlgebraicSet(I::PolynomialIdeal{T, PT}, solver::S) where {T, PT, S} = AlgebraicSet{T, PT, S}(I, solver)
AlgebraicSet{T, PT}() where {T, PT} = AlgebraicSet(PolynomialIdeal{T, PT}(), defaultalgebraicsolver(T))
AlgebraicSet(p::Vector, solver) = AlgebraicSet(PolynomialIdeal(p), solver)

algebraicset(p::Vector, lib::DefaultAlgebraicSetLibrary) = AlgebraicSet(p, lib.solver)
algebraicset(p::Vector, solver) = algebraicset(p, defaultalgebraicsetlibrary(p, solver))
algebraicset(p::Vector) = algebraicset(p, defaultalgebraicsetlibrary(p))

function projectivealgebraicset(p::Vector, lib::DefaultAlgebraicSetLibrary)
    V = AlgebraicSet(p, lib.solver)
    V.projective = true
    V
end
projectivealgebraicset(p::Vector, solver) = projectivealgebraicset(p, defaultalgebraicsetlibrary(p, solver))
projectivealgebraicset(p::Vector) = projectivealgebraicset(p, defaultalgebraicsetlibrary(p))

nequalities(V::AlgebraicSet) = length(V.I.p)
equalities(V::AlgebraicSet) = V.I.p
addequality!(V::AlgebraicSet, p) = push!(V.I.p, p)
Base.intersect(S::AlgebraicSet, T::AlgebraicSet) = AlgebraicSet(S.I + T.I, S.solver)

function Base.show(io::IO, V::AbstractAlgebraicSet)
    print(io, "Algebraic Set defined by ")
    n = nequalities(V)
    if n == 0
        println(io, "no equality")
    elseif n == 1
        println(io, "1 equality")
    else
        println(io, "$n equalities")
    end
    for p in equalities(V)
        println(io, " $p == 0")
    end
end

defaultalgebraicsolver(V::AlgebraicSet) = V.solver
function elements(V::AlgebraicSet{T})::Nullable{Vector{eltype(V)}} where T
    if V.projective
        I = V.I
        els = Vector{T}[]
        for v in variables(I)
            I1 = I + PolynomialIdeal([v - 1])
            els1 = elements(AlgebraicSet(I1, V.solver))
            if isnull(els1)
                return nothing
            end
            append!(els, get(els1))
            I = I + PolynomialIdeal([v])
        end
        els
    else
        solvealgebraicequations(V, V.solver)
    end
end
function computeelements!(V::AlgebraicSet{T}) where T
    if !V.elementscomputed
        els = elements(V)
        V.iszerodimensional = !isnull(els)
        if !isnull(els)
            V.elements = get(els)
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
