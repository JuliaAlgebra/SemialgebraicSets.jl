export iszerodimensional
export algebraicset, projectivealgebraicset, equalities

struct DefaultAlgebraicSetLibrary{S<:AbstractAlgebraicSolver}
    solver::S
end

defaultalgebraicsetlibrary(::Vector{<:APL}, solver::AbstractAlgebraicSolver) = DefaultAlgebraicSetLibrary(solver)
defaultalgebraicsetlibrary(p::Vector{<:APL}, solveroralgo...) = defaultalgebraicsetlibrary(p, defaultalgebraicsolver(p, solveroralgo...))

mutable struct AlgebraicSet{T, PT<:APL{T}, A, S<:AbstractAlgebraicSolver} <: AbstractAlgebraicSet
    I::PolynomialIdeal{T, PT, A}
    projective::Bool
    elements::Vector{Vector{T}}
    elementscomputed::Bool
    iszerodimensional::Bool
    solver::S
end
AlgebraicSet{T, PT, A, S}(I::PolynomialIdeal{T, PT, A}, solver::S) where {T, PT, A, S} = AlgebraicSet{T, PT, A, S}(I, false, Vector{T}[], false, false, solver)
AlgebraicSet(I::PolynomialIdeal{T, PT, A}, solver::S) where {T, PT, A, S} = AlgebraicSet{T, PT, A, S}(I, solver)
AlgebraicSet{T, PT}() where {T, PT} = AlgebraicSet(PolynomialIdeal{T, PT}(), defaultalgebraicsolver(T))
AlgebraicSet(p::Vector, algo::AbstractGröbnerBasisAlgorithm, solver) = AlgebraicSet(ideal(p, algo), solver)

MP.changecoefficienttype(::Type{AlgebraicSet{U, PU, A, S}}, T::Type) where {U, PU, A, S} = AlgebraicSet{T, MP.changecoefficienttype(PU, T), A, S}
function Base.convert(::Type{AlgebraicSet{T, PT, A, S}}, set::AlgebraicSet{T, PT, A, S}) where {T, PT<:APL{T}, A, S<:AbstractAlgebraicSolver}
    return set
end
function Base.convert(::Type{AlgebraicSet{T, PT, A, S}}, set::AlgebraicSet) where {T, PT, A, S}
    return AlgebraicSet{T, PT, A, S}(set.I, set.projective, set.elements, set.elementscomputed, set.iszerodimensional, set.solver)
end

algebraicset(p::Vector, lib::DefaultAlgebraicSetLibrary) = AlgebraicSet(p, defaultgröbnerbasisalgorithm(p), lib.solver)
algebraicset(p::Vector, algo::AbstractGröbnerBasisAlgorithm=defaultgröbnerbasisalgorithm(p), lib::DefaultAlgebraicSetLibrary=defaultalgebraicsetlibrary(p)) = AlgebraicSet(p, algo, lib.solver)
algebraicset(p::Vector, solver) = algebraicset(p, defaultalgebraicsetlibrary(p, solver))

projectivealgebraicset(p::Vector, lib::DefaultAlgebraicSetLibrary) = projectivealgebraicset(p, defaultgröbnerbasisalgorithm(p), lib.solver)
function projectivealgebraicset(p::Vector, algo::AbstractGröbnerBasisAlgorithm=defaultgröbnerbasisalgorithm(p), lib::DefaultAlgebraicSetLibrary=defaultalgebraicsetlibrary(p))
    V = AlgebraicSet(p, algo, lib.solver)
    V.projective = true
    V
end
projectivealgebraicset(p::Vector, algo, solver) = projectivealgebraicset(p, algo, defaultalgebraicsetlibrary(p, solver))
projectivealgebraicset(p::Vector, solver) = projectivealgebraicset(p, defaultgröbnerbasisalgorithm(p), defaultalgebraicsetlibrary(p, solver))

ideal(V::AlgebraicSet) = V.I

MP.variables(V::AlgebraicSet) = MP.variables(V.I)
nequalities(V::AlgebraicSet) = length(V.I.p)
equalities(V::AlgebraicSet) = V.I.p
addequality!(V::AlgebraicSet, p) = push!(V.I.p, p)
Base.intersect(S::AlgebraicSet, T::AlgebraicSet) = AlgebraicSet(S.I + T.I, S.solver)

function algebraicset(set::AbstractAlgebraicSet, solver...)
    return algebraicset(equalities(set), solver...)
end
function Base.show(io::IO, V::AbstractAlgebraicSet)
    print(io, "{ (", join(variables(V), ", "), ") | ",
          join(string.(equalities(V)) .* " = 0", ", "), " }")
end
function _show_els(io::IO, name, n, els, sign)
    if n == 0
        println(io, "no $(name)y")
    elseif n == 1
        println(io, "1 $(name)ty")
    else
        println(io, "$n $(name)ies")
    end
    for p in els
        println(io, " $p $sign 0")
    end
end
function Base.show(io::IO, mime::MIME"text/plain", V::AbstractAlgebraicSet)
    print(io, "Algebraic Set defined by ")
    _show_els(io, "equalit", nequalities(V), equalities(V), "=")
end

defaultalgebraicsolver(V::AlgebraicSet) = V.solver
function elements(V::AlgebraicSet{T}) where T
    if V.projective
        I = V.I
        els = Vector{T}[]
        for v in variables(I)
            I1 = I + ideal([v - 1], V.I.algo)
            els1 = elements(AlgebraicSet(I1, V.solver))
            if els1 === nothing
                return nothing
            end
            append!(els, els1)
            I = I + ideal([v], V.I.algo)
        end
        els
    else
        solvealgebraicequations(V, V.solver)
    end
end
function computeelements!(V::AlgebraicSet{T}) where T
    if !V.elementscomputed
        els = elements(V)
        V.iszerodimensional = els !== nothing
        if V.iszerodimensional
            V.elements = els
        end
        V.elementscomputed = true
    end
end
function iszerodimensional(V::AlgebraicSet)
    computeelements!(V)
    V.iszerodimensional
end

Base.eltype(V::AlgebraicSet{T}) where T = Vector{T}

for f in [:length, :iterate, :lastindex, :getindex]
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
