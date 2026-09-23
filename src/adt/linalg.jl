# 链级线性代数（后端：FiniteMPSAlgorithms；Dense1DTN 上的薄适配层）
#
# dot/norm/distance 的实现委托给 FiniteMPSAlgorithms 在 CanonicalMPS 视图上的
# 同名函数（含 scaling^L 约定，见 fmabackend.jl）。
LinearAlgebra.dot(psiA::Dense1DTN, psiB::Dense1DTN) = dot(psiA.parent, psiB.parent)
LinearAlgebra.norm(psi::Dense1DTN) = norm(psi.parent)

"""
    distance(a::Dense1DTN, b::Dense1DTN)
    distance2(a::Dense1DTN, b::Dense1DTN)

Distance between two tensor networks (`ADT` / `ProcessTensor`), based on the inner products of the site tensors and including the overall `scaling` factors.
`distance = sqrt(distance2)`. Commonly used to verify the accuracy of multiplication/compression.
"""
function distance2(a::Dense1DTN, b::Dense1DTN)
	sA = real(dot(a, a))
	sB = real(dot(b, b))
	c = dot(a, b)
	return abs(sA + sB - 2 * real(c))
end
distance(a::Dense1DTN, b::Dense1DTN) = sqrt(distance2(a, b))


function LinearAlgebra.lmul!(f::Number, psi::Dense1DTN)
    if !isempty(psi)
        psi[1] *= f
    end
    _renormalize!(psi, psi[1], false)
    return psi
end

Base.:*(psi::Dense1DTN, f::Number) = lmul!(f, copy(psi))
Base.:*(f::Number, psi::Dense1DTN) = psi * f
Base.:/(psi::Dense1DTN, f::Number) = psi * (1/f)
Base.:(-)(psi::Dense1DTN) = (-1) * psi


# 精确（未压缩）乘积：物理指标共享的逐点（Hadamard）乘积，
# 委托给 FiniteMPSAlgorithms 的 ⊙（作用在内层 CanonicalMPS 上）
function Base.:*(x::ADT, y::ADT)
    (length(x) == length(y)) || throw(DimensionMismatch())
    @assert !isempty(x)
    return ADT(⊙(x.parent, y.parent))
end

# 块对角直和：FiniteMPSAlgorithms 的实现会把两边的 scaling 折入数据
Base.:+(x::ADT, y::ADT) = ADT(x.parent + y.parent)
Base.:-(x::ADT, y::ADT) = x + (-y)


permute!(x::ADT, perm::AbstractVector{Int}; kwargs...) = (permute!(x.parent, perm; kwargs...); x)
permute(x::ADT, perm::AbstractVector{Int}; kwargs...) = ADT(permute(x.parent, perm; kwargs...))

function _mult_site_n(xj::DenseMPSTensor, yj::DenseMPSTensor)
    @tensor r[1,4,2,5;3,6] := xj[1,2,3] * yj[4,5,6]
    return r
end
