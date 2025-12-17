abstract type HybridizationStyle end


struct AdditiveBath{T<:Number} <: HybridizationStyle
	op::Vector{T}
end

phydim(b::AdditiveBath) = length(b.op)

struct NonAdditiveBath{T<:Number} <: HybridizationStyle
	op::Matrix{T}

function NonAdditiveBath{T}(op::AbstractMatrix) where {T<:Number}
	(size(op, 1) == size(op, 2)) || throw(ArgumentError("square matrix expected"))
	new{T}(convert(Matrix{T}, op))
end
end
NonAdditiveBath(a::AbstractMatrix{T}) where {T<:Number} = NonAdditiveBath{T}(a)

phydim(b::NonAdditiveBath) = size(b.op, 1)

include("partialif/partialif.jl")