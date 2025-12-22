
"""
	ProcessTensor{A <: MPOTensor}
Finite Matrix Product Operator which stores a chain of rank-4 site tensors.
"""
struct ProcessTensor{T<:Number} <: Dense1DTN{T}
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
function ProcessTensor{T}(mpotensors::AbstractVector, scaling::Ref{Float64}) where {T<:Number}
	_check_mpo_space(mpotensors)
	return new{T}(convert(Vector{Array{T, 4}}, mpotensors), scaling)
end

end
ProcessTensor(data::AbstractVector{<:DenseMPOTensor{T}}; scaling::Real=1) where {T <: Number} = ProcessTensor{T}(data, Ref(float(scaling)))

function ProcessTensor(::Type{T}, ds::AbstractVector{Int}) where {T <: Number}
	data = [reshape(_eye(T, d), 1, d, 1, d) for d in ds]
	return ProcessTensor(data, scaling=1)
end
ProcessTensor(ds::AbstractVector{Int}) = ProcessTensor(Float64, ds)
ProcessTensor(::Type{T}, L::Int; d::Int=2) where {T <: Number} = ProcessTensor(T, [d for i in 1:L])
ProcessTensor(L::Int; d::Int=2) = ProcessTensor(Float64, L, d=d)

Base.copy(psi::ProcessTensor) = ProcessTensor(copy(psi.data), scaling=scaling(psi))


isleftcanonical(a::ProcessTensor; kwargs...) = all(x->isleftcanonical(x; kwargs...), a.data)
isrightcanonical(a::ProcessTensor; kwargs...) = all(x->isrightcanonical(x; kwargs...), a.data)


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
	return true
end


# initializers

"""
	randompt(::Type{T}, ds::Vector{Int}; D::Int) where {T<:Number}
	dy are the input dimensions, dx are the output dimensions
"""
function randompt(::Type{T}, ds::Vector{Int}; D::Int) where {T<:Number}
	L = length(ds)
	r = Vector{Array{T, 4}}(undef, L)
	r[1] = randn(T, 1, ds[1], D, ds[1])
	r[L] = randn(T, D, ds[L], 1, ds[L])
	for i in 2:L-1
		r[i] = randn(T, D, ds[i], D, ds[i])
	end
	return ProcessTensor(r)
end 
randompt(::Type{T}, L::Int; d::Int=2, D::Int) where {T<:Number} = randompt(T, [d for i in 1:L], D=D)
randompt(L::Int; kwargs...) = randompt(Float64, L; kwargs...)



