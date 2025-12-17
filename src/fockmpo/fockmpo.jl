
"""
	FockMPO{A <: MPOTensor}
Finite Matrix Product Operator which stores a chain of rank-4 site tensors.
"""
struct FockMPO{T<:Number} <: Dense1DTN{T}
	data::Vector{Array{T, 4}}
	scaling::Ref{Float64}
"""
	DenseMPO{A}(mpotensors::Vector)
Constructor entrance for MPO, which only supports strictly quantum number conserving operators

site tensor convention:
i mean in arrow, o means out arrow
    o 
    |
    2
o-1   3-i
	4
	|
	i
The left and right boundaries are always vacuum.
The case that the right boundary is not vacuum corresponds to operators which do not conserve quantum number, 
such as a†, this case is implemented with another MPO object.
"""
function FockMPO{T}(mpotensors::AbstractVector, scaling::Ref{Float64}) where {T<:Number}
	_check_mpo_space(mpotensors)
	return new{T}(convert(Vector{Array{T, 4}}, mpotensors), scaling)
end

end
FockMPO(data::AbstractVector{<:DenseMPOTensor{T}}; scaling::Real=1) where {T <: Number} = MPO{T}(data, Ref(float(scaling)))


# attributes

function _check_mpo_space(mpotensors::Vector)
	L = length(mpotensors)
	for i in 1:L-1
		(space_r(mpotensors[i]) == space_l(mpotensors[i+1])) || throw(DimensionMismatch())
	end
	for m in mpotensors
		(size(m, 2) == size(m, 4)) || throw(ArgumentError("physical dimension mismatch"))
	end
	# boundaries should be dimension 
	(space_l(mpotensors[1]) == 1) || throw(DimensionMismatch())
	(space_r(mpotensors[L]) == 1) || throw(DimensionMismatch())
	all(x->size(x, 2)==size(x,4)==2, mpotensors) || throw(ArgumentError("physical dimension of site tensors must be 2"))
	return true
end


# initializers

"""
	randomfockmpo(::Type{T}, ds::Vector{Int}; D::Int) where {T<:Number}
	dy are the input dimensions, dx are the output dimensions
"""
function randomfockmpo(::Type{T}, ds::Vector{Int}; D::Int) where {T<:Number}
	L = length(dx)
	r = Vector{Array{T, 4}}(undef, L)
	r[1] = randn(T, 1, ds[1], D, ds[1])
	r[L] = randn(T, D, ds[L], 1, ds[L])
	for i in 2:L-1
		r[i] = randn(T, D, ds[i], D, ds[i])
	end
	return FockMPO(r)
end 
randomfockmpo(::Type{T}, L::Int; d::Int, D::Int) where {T<:Number} = randomfockmpo(T, [d for i in 1:L], D=D)
randomfockmpo(L::Int; kwargs...) = randomfockmpo(Float64, L; kwargs...)


include("linalg.jl")
include("orth.jl")
include("integrate.jl")
