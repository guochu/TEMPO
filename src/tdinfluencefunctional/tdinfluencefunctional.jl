
abstract type TdHybridizationStyle <: HybridizationStyle end
(x::TdHybridizationStyle)(t::Real) = x.op(t)

"""
    struct AdditiveTdHyb{F}

    op should be a vector function, 
    i.e., op(t) gives a vector for any t
"""
struct AdditiveTdHyb{F<:Function} <: TdHybridizationStyle
    op::F
end
AdditiveTdHyb(x::AbstractVector{<:Real}, f::Function) = AdditiveTdHyb(t->x .* f(t))

phydim(b::AdditiveTdHyb) = length(b(0))

TO.scalartype(::Type{AdditiveTdHyb}) = Float64
function pairop(b::AdditiveTdHyb, t::Real)
	op = Matrix(Diagonal(b(t)))
	return op, op
end
function AdditiveTdHyb(a::AbstractMatrix, f::Function)
	isdiag(a) || throw(ArgumentError("AdditiveHyb requires diagonal matrix"))
	adiag = [a[i, i] for i in 1:size(a, 1)]
	return AdditiveTdHyb(adiag, f)
end

abstract type GeneralTdHybStyle <: TdHybridizationStyle end

"""
	struct NonAdditiveTdHyb{F<:Function}

The Impurity-bath coupling is op*(a+a')
"""
struct NonAdditiveTdHyb{F<:Function} <: GeneralTdHybStyle
	op::F
end
function NonAdditiveTdHyb(op::AbstractMatrix, f::Function) 
	(size(op, 1) == size(op, 2)) || throw(ArgumentError("square matrix expected"))
	ishermitian(op) || throw(ArgumentError("Hermitian matrix required"))
	return NonAdditiveTdHyb(t->op .* f(t))
end
phydim(b::NonAdditiveHyb) = size(b(0), 1)

TO.scalartype(t::NonAdditiveTdHyb) = scalartype(t(0.))
function pairop(b::NonAdditiveTdHyb, t::Real)
    r = b(t)
    return r, r
end 

"""
	struct NonDiagonalTdHyb{T<:Number}

The Impurity-bath coupling is op*a + op'*a'
"""
struct NonDiagonalTdHyb{F<:Function} <: GeneralTdHybStyle
	op::F
end
function NonDiagonalTdHyb(op::AbstractMatrix, f::Function) 
	(size(op, 1) == size(op, 2)) || throw(ArgumentError("square matrix expected"))
	return NonAdditiveTdHyb(t->op .* f(t))
end
NonDiagonalTdHyb(hyb::NonAdditiveTdHyb) = NonDiagonalTdHyb(hyb.op)

phydim(b::NonDiagonalTdHyb) = size(b(0), 1)

TO.scalartype(t::NonAdditiveTdHyb) = scalartype(t(0.))
function pairop(b::NonAdditiveTdHyb, t::Real)
    r = b(t)
    return r, r'
end 

include("partialif/partialif.jl")
include("ttiif/ttiif.jl")