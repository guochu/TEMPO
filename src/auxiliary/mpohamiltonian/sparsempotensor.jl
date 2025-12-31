"""
	struct SparseMPOTensor{M<:MPOTensor, T<:Number} <: AbstractSparseMPOTensor{S}

leftspaces are left space, rightspaces are right space
"""
struct SparseMPOTensor{M<:AbstractMatrix, T<:Number} <: AbstractSparseMPOTensor
	Os::Array{Union{M, T}, 2}
	d::Int
end
Base.copy(x::SparseMPOTensor) = SparseMPOTensor(copy(x.Os), phydim(x))
TO.scalartype(::Type{SparseMPOTensor{M,T}}) where {M,T} = T

function SparseMPOTensor{M, T}(data::AbstractMatrix) where {M <:AbstractMatrix, T<:Number}
	Os, pspace = compute_mpotensor_data(M, T, data)
	return SparseMPOTensor{M, T}(Os, pspace)
end


"""
	SparseMPOTensor(data::Array{Any, 2}) 
"""
function SparseMPOTensor(data::AbstractMatrix) 
	T = compute_eltype(data)
	M = Matrix{T}
	return SparseMPOTensor{M, T}(data)
end
SparseMPOTensor(data::AbstractMatrix{Union{M, T}}) where {M<:AbstractMatrix, T<:Number} = SparseMPOTensor{M, T}(data)


Base.contains(m::SparseMPOTensor{M, T}, i::Int, j::Int) where {M, T} = (m.Os[i, j] != zero(T))
function isscal(x::SparseMPOTensor{M,T}, i::Int, j::Int) where {M,T}
	sj = x.Os[i, j]
	return (sj isa T) && (abs(sj) > 1.0e-14)
end 

Base.getindex(m::SparseMPOTensor, i::Union{UnitRange, Colon}, j::Union{UnitRange, Colon}) = SparseMPOTensor(m.Os[i, j], phydim(m))
Base.getindex(m::SparseMPOTensor, i::Int, j::Union{UnitRange, Colon}) = getindex(m, i:i, j)
Base.getindex(m::SparseMPOTensor, i::Union{UnitRange, Colon}, j::Int) = getindex(m, i, j:j)

function Base.setindex!(m::SparseMPOTensor{M, T}, v, i::Int, j::Int) where {M, T}
	if isa(v, Number)
		m.Os[i, j] = convert(T, v)
	elseif isa(v, MPSBondTensor)
		(size(v, 1) == size(v, 2) == phydim(m)) || throw(DimensionMismatch("input matrix size mismatch with phydim"))
		m.Os[i, j] = convert(M, v)
	else
		throw(ArgumentError("input should be scalar or Matrix type"))
	end
end