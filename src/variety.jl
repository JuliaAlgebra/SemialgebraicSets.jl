export is_zero_dimensional
export algebraic_set, projective_algebraic_set, equalities

struct DefaultAlgebraicSetLibrary{S<:AbstractAlgebraicSolver}
    solver::S
end

function default_gröbner_basis_algorithm(p, lib::DefaultAlgebraicSetLibrary)
    return default_gröbner_basis_algorithm(p, lib.solver)
end

function default_algebraic_set_library(
    ::Vector{<:_APL},
    solver::AbstractAlgebraicSolver,
)
    return DefaultAlgebraicSetLibrary(solver)
end
function default_algebraic_set_library(p::Vector{<:_APL}, solveroralgo...)
    return default_algebraic_set_library(
        p,
        default_algebraic_solver(p, solveroralgo...),
    )
end

mutable struct AlgebraicSet{T,PT<:_APL{T},A,S<:AbstractAlgebraicSolver,U} <:
               AbstractAlgebraicSet
    I::PolynomialIdeal{T,PT,A}
    projective::Bool
    elements::Vector{Vector{U}}
    elements_computed::Bool
    is_zero_dimensional::Bool
    solver::S
end
function AlgebraicSet{T,PT,A,S,U}(
    I::PolynomialIdeal{T,PT,A},
    solver::S,
) where {T,PT,A,S,U}
    return AlgebraicSet{T,PT,A,S,U}(I, false, Vector{U}[], false, false, solver)
end
function AlgebraicSet(I::PolynomialIdeal{T,PT,A}, solver::S) where {T,PT,A,S}
    return AlgebraicSet{T,PT,A,S,promote_for(T, S)}(I, solver)
end
function AlgebraicSet{T,PT}() where {T,PT}
    return AlgebraicSet(PolynomialIdeal{T,PT}(), default_algebraic_solver(T))
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
) where {T,PT<:_APL{T},A,S<:AbstractAlgebraicSolver,U}
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
        set.elements_computed,
        set.is_zero_dimensional,
        set.solver,
    )
end

function algebraic_set(
    p::Vector,
    algo::AbstractGröbnerBasisAlgorithm = default_gröbner_basis_algorithm(p),
    lib::DefaultAlgebraicSetLibrary = default_algebraic_set_library(p, algo);
    projective::Bool = false,
)
    set = AlgebraicSet(p, algo, lib.solver)
    set.projective = projective
    return set
end

function algebraic_set(p::Vector, lib::DefaultAlgebraicSetLibrary; kws...)
    return algebraic_set(
        p,
        default_gröbner_basis_algorithm(p, lib),
        lib.solver;
        kws...,
    )
end

function algebraic_set(p::Vector, solver; kws...)
    return algebraic_set(p, default_algebraic_set_library(p, solver); kws...)
end

function algebraic_set(
    p::Vector,
    algo::AbstractGröbnerBasisAlgorithm,
    solver;
    kws...,
)
    return algebraic_set(
        p,
        algo,
        default_algebraic_set_library(p, solver);
        kws...,
    )
end

function projective_algebraic_set(p::Vector, args...)
    return algebraic_set(p, args...; projective = true)
end

ideal(V::AlgebraicSet) = V.I
function MA.promote_operation(
    ::typeof(ideal),
    ::Type{<:AlgebraicSet{T,P,A}},
) where {T,P,A}
    return PolynomialIdeal{T,P,A}
end

MP.variables(V::AlgebraicSet) = MP.variables(V.I)
MP.monomial_type(::Type{<:AlgebraicSet{T,P}}) where {T,P} = MP.monomial_type(P)
nequalities(V::AlgebraicSet) = length(V.I.p)
equalities(V::AlgebraicSet) = V.I.p
add_equality!(V::AlgebraicSet, p) = push!(V.I.p, p)
function Base.intersect(S::AlgebraicSet, T::AlgebraicSet)
    return AlgebraicSet(S.I + T.I, S.solver)
end

function algebraic_set(set::AbstractAlgebraicSet, solver...)
    return algebraic_set(equalities(set), solver...)
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

default_algebraic_solver(V::AlgebraicSet) = V.solver
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
        return els
    else
        return solve(V, V.solver)
    end
end
function compute_elements!(V::AlgebraicSet{T}) where {T}
    if !V.elements_computed
        els = elements(V)
        V.is_zero_dimensional = els !== nothing
        if V.is_zero_dimensional
            V.elements = els
        end
        V.elements_computed = true
    end
end
function is_zero_dimensional(V::AlgebraicSet)
    compute_elements!(V)
    return V.is_zero_dimensional
end

Base.eltype(V::AlgebraicSet{T,PT,A,S,U}) where {T,PT,A,S,U} = Vector{U}

for f in [:length, :iterate, :lastindex, :getindex]
    @eval begin
        function Base.$f(V::AlgebraicSet, args...)
            compute_elements!(V)
            if !is_zero_dimensional(V)
                error("A non zero-dimensional algebraic set is not iterable")
            end
            return $f(V.elements, args...)
        end
    end
end
