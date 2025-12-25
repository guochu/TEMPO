abstract type HybridizationStyle end


struct AdditiveHyb <: HybridizationStyle
	op::Vector{Float64}
end

AdditiveHyb(x::AbstractVector{<:Real}) = AdditiveHyb(float(x))

phydim(b::AdditiveHyb) = length(b.op)

TO.scalartype(::Type{AdditiveHyb}) = Float64

struct NonAdditiveHyb{T<:Number} <: HybridizationStyle
	op::Matrix{T}

function NonAdditiveHyb{T}(op::AbstractMatrix) where {T<:Number}
	(size(op, 1) == size(op, 2)) || throw(ArgumentError("square matrix expected"))
	new{T}(convert(Matrix{T}, op))
end
end
NonAdditiveHyb(a::AbstractMatrix{T}) where {T<:Number} = NonAdditiveHyb{T}(a)

phydim(b::NonAdditiveHyb) = size(b.op, 1)

TO.scalartype(::Type{NonAdditiveHyb{T}}) where T = T


# struct NonDiagonalHyb{T<:Number} <: HybridizationStyle
# 	sp::Matrix{T}
# end

# phydim(b::NonDiagonalHyb) = size(b.op, 1)

# TO.scalartype(::Type{NonDiagonalHyb{T}}) where T = T


include("partialif/partialif.jl")
include("ptpartialif/ptpartialif.jl")