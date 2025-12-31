

# coeff * α^n, α must be in [0, 1]
"""
	struct ExponentialDecayTerm{M1, M, M2, T <:Number}

Exponential decay of the form coeff * [â ⊗ (αm̂)^⊗n ⊗ b̂]
"""
struct ExponentialDecayTerm{M1<:AbstractMatrix, M<:AbstractMatrix, M2, T <:Number} <: AbstractLongRangeTerm
    a::M1
    m::M
    b::M2
    α::T
    coeff::T
end

function ExponentialDecayTerm(a::AbstractMatrix, b::AbstractMatrix; middle::AbstractMatrix=_eye(size(a, 1)), α::Number=1., coeff::Number=1.) 
    T = promote_type(typeof(α), typeof(coeff))
    return ExponentialDecayTerm(a, middle, b, convert(T, α), convert(T, coeff))
end

TO.scalartype(::Type{ExponentialDecayTerm{M1, M, M2, T}}) where {M1, M, M2, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), T)

Base.adjoint(x::ExponentialDecayTerm) = ExponentialDecayTerm(_op_adjoint(x.a, x.m, x.b)..., conj(x.α), conj(coeff(x)))
_op_adjoint(a::AbstractMatrix, m::AbstractMatrix, b::AbstractMatrix) = (a', m', b')


function _longrange_schurmpo_util(h1, h2s::Vector{<:ExponentialDecayTerm})
    isempty(h2s) && throw(ArgumentError("empty interactions."))
	pspace = size(h2s[1].a, 1)
	N = length(h2s)
	T = Float64
	for item in h2s
		T = promote_type(T, scalartype(item))
	end
	cell = Matrix{Any}(undef, N+2, N+2)
	for i in 1:length(cell)
		cell[i] = zero(T)
	end
	# diagonals
	cell[1, 1] = 1
	cell[end, end] = 1
	cell[1, end] = h1
	for i in 1:N
		cell[i+1, i+1] = h2s[i].α * h2s[i].m
		cell[1, i+1] = h2s[i].coeff * h2s[i].a
		cell[i+1, end] = h2s[i].α * h2s[i].b
	end
	return SchurMPOTensor(cell)
end



"""
    SchurMPOTensor(h1::ScalarSiteOp, h2s::Vector{<:ExponentialDecayTerm})
    SchurMPOTensor(h2s::Vector{<:ExponentialDecayTerm})

Return an SchurMPOTensor, with outer matrix size (N+2)×(N+2) (N=length(h2s))
Algorithm reference: "Time-evolving a matrix product state with long-ranged interactions"
"""
SchurMPOTensor(h1::AbstractMatrix{<:Number}, h2s::Vector{<:ExponentialDecayTerm}) = _longrange_schurmpo_util(h1, h2s)
SchurMPOTensor(h2s::Vector{<:ExponentialDecayTerm}) = _longrange_schurmpo_util(0., h2s)
