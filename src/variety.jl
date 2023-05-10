export iszerodimensional
export algebraicset, projectivealgebraicset, equalities

struct DefaultAlgebraicSetLibrary{S<:AbstractAlgebraicSolver}
    solver::S
end

function defaultalgebraicsetlibrary(
    ::Vector{<:APL},
    solver::AbstractAlgebraicSolver,
)
    return DefaultAlgebraicSetLibrary(solver)
end
function defaultalgebraicsetlibrary(p::Vector{<:APL}, solveroralgo...)
    return defaultalgebraicsetlibrary(
        p,
        defaultalgebraicsolver(p, solveroralgo...),
    )
end

mutable struct AlgebraicSet{T,PT<:APL{T},A,S<:AbstractAlgebraicSolver,U} <:
               AbstractAlgebraicSet
    I::PolynomialIdeal{T,PT,A}
    projective::Bool
    elements::Vector{Vector{U}}
    elementscomputed::Bool
    iszerodimensional::Bool
    solver::S
end
function AlgebraicSet{T,PT,A,S,U}(
    I::PolynomialIdeal{T,PT,A},
    solver::S,
) where {T,PT,A,S,U}
    return AlgebraicSet{T,PT,A,S,U}(I, false, Vector{U}[], false, false, solver)
end
function AlgebraicSet(I::PolynomialIdeal{T,PT,A}, solver::S) where {T,PT,A,S}
    return AlgebraicSet{T,PT,A,S,float(T)}(I, solver)
end
function AlgebraicSet{T,PT}() where {T,PT}
    return AlgebraicSet(PolynomialIdeal{T,PT}(), defaultalgebraicsolver(T))
end
function AlgebraicSet(p::Vector, algo::AbstractGröbnerBasisAlgorithm, solver)
    return AlgebraicSet(ideal(p, algo), solver)
end

function MP.similar_type(
    ::Type{AlgebraicSet{U,PU,A,S,UU}},
    T::Type,
) where {U,PU,A,S,UU}
    return AlgebraicSet{T,MP.similar_type(PU, T),A,S,float(T)}
end
function Base.convert(
    ::Type{AlgebraicSet{T,PT,A,S,U}},
    set::AlgebraicSet{T,PT,A,S,U},
) where {T,PT<:APL{T},A,S<:AbstractAlgebraicSolver,U}
    return set
end
function Base.convert(
    ::Type{AlgebraicSet{T,PT,A,S,U}},
    set::AlgebraicSet,
) where {T,PT,A,S,U}
    return AlgebraicSet{T,PT,A,S,U}(
        set.I,
        set.projective,
        set.elements,
        set.elementscomputed,
        set.iszerodimensional,
        set.solver,
    )
end

function algebraicset(p::Vector, lib::DefaultAlgebraicSetLibrary)
    return AlgebraicSet(p, defaultgröbnerbasisalgorithm(p), lib.solver)
end
function algebraicset(
    p::Vector,
    algo::AbstractGröbnerBasisAlgorithm = defaultgröbnerbasisalgorithm(p),
    lib::DefaultAlgebraicSetLibrary = defaultalgebraicsetlibrary(p),
)
    return AlgebraicSet(p, algo, lib.solver)
end
function algebraicset(p::Vector, solver)
    return algebraicset(p, defaultalgebraicsetlibrary(p, solver))
end
function algebraicset(p::Vector, algo::AbstractGröbnerBasisAlgorithm, solver)
    return algebraicset(p, algo, defaultalgebraicsetlibrary(p, solver))
end

function projectivealgebraicset(p::Vector, lib::DefaultAlgebraicSetLibrary)
    return projectivealgebraicset(
        p,
        defaultgröbnerbasisalgorithm(p),
        lib.solver,
    )
end
function projectivealgebraicset(
    p::Vector,
    algo::AbstractGröbnerBasisAlgorithm = defaultgröbnerbasisalgorithm(p),
    lib::DefaultAlgebraicSetLibrary = defaultalgebraicsetlibrary(p),
)
    V = AlgebraicSet(p, algo, lib.solver)
    V.projective = true
    return V
end
function projectivealgebraicset(p::Vector, algo, solver)
    return projectivealgebraicset(
        p,
        algo,
        defaultalgebraicsetlibrary(p, solver),
    )
end
function projectivealgebraicset(p::Vector, solver)
    return projectivealgebraicset(
        p,
        defaultgröbnerbasisalgorithm(p),
        defaultalgebraicsetlibrary(p, solver),
    )
end

ideal(V::AlgebraicSet) = V.I

MP.variables(V::AlgebraicSet) = MP.variables(V.I)
nequalities(V::AlgebraicSet) = length(V.I.p)
equalities(V::AlgebraicSet) = V.I.p
addequality!(V::AlgebraicSet, p) = push!(V.I.p, p)
function Base.intersect(S::AlgebraicSet, T::AlgebraicSet)
    return AlgebraicSet(S.I + T.I, S.solver)
end

function algebraicset(set::AbstractAlgebraicSet, solver...)
    return algebraicset(equalities(set), solver...)
end
function Base.show(io::IO, V::AbstractAlgebraicSet)
    return print(
        io,
        "{ (",
        join(variables(V), ", "),
        ") | ",
        join(string.(equalities(V)) .* " = 0", ", "),
        " }",
    )
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
    return _show_els(io, "equalit", nequalities(V), equalities(V), "=")
end

defaultalgebraicsolver(V::AlgebraicSet) = V.solver
function elements(V::AlgebraicSet{T,PT,A,S,U}) where {T,PT,A,S,U}
    if V.projective
        I = V.I
        els = Vector{U}[]
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
function computeelements!(V::AlgebraicSet{T}) where {T}
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
    return V.iszerodimensional
end

Base.eltype(V::AlgebraicSet{T,PT,A,S,U}) where {T,PT,A,S,U} = Vector{U}

for f in [:length, :iterate, :lastindex, :getindex]
    @eval begin
        function Base.$f(V::AlgebraicSet, args...)
            computeelements!(V)
            if !iszerodimensional(V)
                error("A non zero-dimensional algebraic set is not iterable")
            end
            return $f(V.elements, args...)
        end
    end
end
