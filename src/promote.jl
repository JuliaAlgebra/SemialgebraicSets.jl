function SA.promote_bases_with_maps(f::FullSpace, p::Union{MP.AbstractPolynomialLike,AbstractSemialgebraicSet})
    return (f, nothing), (p, nothing)
end

function SA.promote_bases_with_maps(a::AbstractSemialgebraicSet, b::Union{MP.AbstractPolynomialLike,AbstractSemialgebraicSet})
    _a, _b = MP.promote_variables_with_maps(MP.variables(a), MP.variables(b))
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end
