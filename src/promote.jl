function _promote_variables_with_maps(vars_a, vars_b)
    if vars_a == vars_b
        return (vars_a, nothing), (vars_b, nothing)
    end
    all_vars = sort!(union(vars_a, vars_b))
    return (all_vars, all_vars), (all_vars, all_vars)
end

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
    _a, _b = _promote_variables_with_maps(MP.variables(a), MP.variables(b))
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end

function _promote_poly(p, new_vars)
    return p * prod(v -> v^0, new_vars)
end

function _promote_polys(polys, new_vars)
    return [_promote_poly(p, new_vars) for p in polys]
end

function SA.promote_with_map(
    p::MP.AbstractPolynomialLike,
    new_vars,
    exponent_map,
)
    return _promote_poly(p, new_vars), exponent_map
end

function SA.promote_with_map(set::AlgebraicSet, new_vars, exponent_map)
    new_polys = _promote_polys(equalities(set), new_vars)
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
    new_p = _promote_polys(inequalities(set), new_vars)
    return BasicSemialgebraicSet(new_V, new_p), exponent_map
end
