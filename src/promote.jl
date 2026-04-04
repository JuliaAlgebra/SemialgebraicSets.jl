const _PolynomialLike = Union{MP.AbstractPolynomialLike,SA.AlgebraElement}

function SA.promote_bases_with_maps(
    ::FullSpace,
    p::Union{MP.AbstractPolynomialLike,AbstractSemialgebraicSet},
)
    return (FullSpace(), nothing), (p, nothing)
end

function SA.promote_bases_with_maps(::FullSpace, p::SA.AlgebraElement)
    return (FullSpace(), nothing), (p, nothing)
end

function SA.promote_bases_with_maps(
    a::AbstractSemialgebraicSet,
    b::Union{MP.AbstractPolynomialLike,AbstractSemialgebraicSet},
)
    _a, _b = MP.promote_variables_with_maps(MP.variables(a), MP.variables(b))
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end

function _promote_ae(b::SA.AlgebraElement, all_vars)
    alg = SA.parent(b)
    new_obj = typeof(SA.object(alg))(all_vars)
    new_basis, _ = SA.promote_with_map(SA.basis(alg), new_obj, identity)
    (new_alg, alg_map), _ = SA.promote_bases_with_maps(alg, new_basis)
    return SA.maybe_promote(b, new_alg, alg_map)
end

function SA.promote_bases_with_maps(
    a::AbstractSemialgebraicSet,
    b::SA.AlgebraElement,
)
    _a, _b = MP.promote_variables_with_maps(MP.variables(a), MP.variables(b))
    return SA.maybe_promote(a, _a...), _promote_ae(b, _b[1])
end

function _map_polys(p, vars, map)
    return first.(SA.promote_with_map.(p, Ref(vars), Ref(map)))
end

function SA.promote_with_map(set::AlgebraicSet, new_vars, exponent_map)
    new_polys = _map_polys(equalities(set), new_vars, exponent_map)
    new_I = PolynomialIdeal(new_polys, set.I.algo)
    new_set = AlgebraicSet(new_I, set.solver)
    return new_set, exponent_map
end

function SA.promote_with_map(::FullSpace, _vars, exponent_map)
    return FullSpace(), exponent_map
end

function SA.promote_with_map(set::FixedVariablesSet, _vars, exponent_map)
    return set, exponent_map
end

function SA.promote_with_map(set::BasicSemialgebraicSet, new_vars, exponent_map)
    new_V, _ = SA.promote_with_map(set.V, new_vars, exponent_map)
    new_polys = _map_polys(inequalities(set), new_vars, exponent_map)
    return BasicSemialgebraicSet(new_V, new_polys), exponent_map
end
