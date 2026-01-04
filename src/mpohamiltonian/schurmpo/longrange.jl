abstract type AbstractLongRangeTerm end

space_l(x::AbstractLongRangeTerm) = isa(x.a, AbstractMatrix) ? size(x.a, 1) : 1
space_r(x::AbstractLongRangeTerm) = isa(x.b, AbstractMatrix) ? size(x.b, 2) : 1
coeff(x::AbstractLongRangeTerm) = x.coeff
