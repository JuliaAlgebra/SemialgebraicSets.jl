function SA.promote_bases_with_maps(
    ::FullSpace,
    p::Union{MP.AbstractPolynomialLike,AbstractSemialgebraicSet},
)
    return (FullSpace(), nothing), (p, nothing)
end

function SA.promote_bases_with_maps(
    a::AbstractSemialgebraicSet,
    b::Union{MP.AbstractPolynomialLike,AbstractSemialgebraicSet},
)
    _a, _b = MP.promote_variables_with_maps(MP.variables(a), MP.variables(b))
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end

function _map_polys(p, vars, map)
    return first.(SA.promote_with_map.(p, Ref(vars), Ref(map)))
end

function _promote_polys_to_common_variables(polys)
    isempty(polys) && return polys
    all_vars = MP.variables(polys)
    return [first(SA.maybe_promote(p, MP.promote_variables_with_maps(MP.variables(p), all_vars)[1]...)) for p in polys]
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
