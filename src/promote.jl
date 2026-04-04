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

function _promote_polys(polys, vars, exponent_map)
    return [_promote_poly(p, vars, exponent_map) for p in polys]
end

function _promote_poly(p::MP.AbstractPolynomialLike, vars, exponent_map)
    return sum(
        MP.coefficient(t) * first(SA.promote_with_map(MP.monomial(t), vars, exponent_map))
        for t in MP.terms(p)
    )
end

function SA.promote_with_map(set::AlgebraicSet, vars, exponent_map)
    new_polys = _promote_polys(equalities(set), vars, exponent_map)
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

function SA.promote_with_map(set::BasicSemialgebraicSet, vars, exponent_map)
    new_V, _ = SA.promote_with_map(set.V, vars, exponent_map)
    new_p = _promote_polys(inequalities(set), vars, exponent_map)
    return BasicSemialgebraicSet(new_V, new_p), exponent_map
end
