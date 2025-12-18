abstract type HybridizationStyle end


struct AdditiveHyb <: HybridizationStyle
	op::Vector{Float64}
end

AdditiveHyb(x::AbstractVector{<:Real}) = AdditiveHyb(float(x))

phydim(b::AdditiveHyb) = length(b.op)

struct NonAdditiveHyb{T<:Number} <: HybridizationStyle
	op::Matrix{T}

function NonAdditiveHyb{T}(op::AbstractMatrix) where {T<:Number}
	(size(op, 1) == size(op, 2)) || throw(ArgumentError("square matrix expected"))
	new{T}(convert(Matrix{T}, op))
end
end
NonAdditiveHyb(a::AbstractMatrix{T}) where {T<:Number} = NonAdditiveHyb{T}(a)

phydim(b::NonAdditiveHyb) = size(b.op, 1)

include("partialif/partialif.jl")