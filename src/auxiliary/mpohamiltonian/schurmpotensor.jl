"""
	SchurMPOTensor{M<:MPOTensor, T<:Number} 

Upper triangular matrix, with each entry an MPOTensor
The top left and bottom right must be identity
"""
struct SchurMPOTensor{M<:AbstractMatrix, T<:Number} <: AbstractSparseMPOTensor
	Os::Array{Union{M, T}, 2}
	d::Int
end

function SchurMPOTensor{M, T}(data::AbstractMatrix) where {M <:AbstractMatrix, T<:Number}
	(size(data, 1) == size(data, 2)) || throw(ArgumentError("SchurMPOTensor requires a square matrix"))
	Os, pspace = compute_mpotensor_data(M, T, data)
	for i in 1:size(Os, 1)
		for j in 1:i-1
			(Os[i, j] == zero(T)) || throw(ArgumentError("SchurMPOTensor should be upper triangular"))
		end
	end
	return SchurMPOTensor{M, T}(Os, pspace)
end

TO.scalartype(::Type{SchurMPOTensor{M,T}}) where {M,T} = T
Base.copy(x::SchurMPOTensor) = SchurMPOTensor(copy(x.Os), phydim(x))

# upper triangular form
# the middle diagonal terms may be identity operator or the JW operator,
"""
	SchurMPOTensor(data::Array{Any, 2}) 
"""
function SchurMPOTensor(data::AbstractMatrix{Any}) 
	T = compute_eltype(data)
	M = Matrix{T}
	return SchurMPOTensor{M, T}(data)
end
SchurMPOTensor(data::AbstractMatrix{SparseMPOTensorElement{M, T}}) where {M<:AbstractMatrix, T<:Number} = SchurMPOTensor{spacetype(M), M, T}(data)


Base.contains(m::SchurMPOTensor{M, T}, i::Int, j::Int) where {M, T} = (i <= j) && (m.Os[i, j] != zero(T))
function isscal(x::SchurMPOTensor{M,T}, i::Int, j::Int) where {M,T}
	sj = x.Os[i, j]
	return (sj isa T) && (abs(sj) > 1.0e-14)
end 

Base.convert(::Type{<:SparseMPOTensor}, t::SchurMPOTensor) = SparseMPOTensor(t.Os, phydim(t))


function Base.setindex!(m::SchurMPOTensor{M, T}, v, i::Int, j::Int) where {M, T}
	(i > j) && throw(ArgumentError("not allowed to set the low triangular portion"))
	if isa(v, Number)
		m.Os[i, j] = convert(T, v)
	elseif isa(v, AbstractMatrix)
		(size(v, 1) == size(v, 2) == phydim(m)) || throw(DimensionMismatch("input matrix size mismatch with phydim"))
		m.Os[i, j] = convert(M, v)
	else
		throw(ArgumentError("input should be scalar or Matrix type"))
	end
end

function compute_eltype(data::AbstractMatrix)
	T = Float64
	for sj in data
		if isa(sj, Number)
			T = promote_type(T, typeof(sj))
		elseif isa(sj, AbstractMatrix)
			T = promote_type(T, scalartype(sj))
		else
			throw(ArgumentError("eltype must be scalar or Matrix"))
		end		
	end
	return T	
end
# the same row should have the same left space
# the same column should have the same right space
function compute_mpotensor_data(::Type{M}, ::Type{T}, data::AbstractMatrix) where {M<:AbstractMatrix, T<:Number}
	@assert !isempty(data)
	m, n = size(data)
	new_data = Array{Union{M, T}, 2}(undef, m, n) 

	# check spaces
	d::Union{Nothing, Int} = nothing
	for i in 1:m
		for j in 1:n
			sj = data[i, j]
 			if isa(sj, AbstractMatrix)
 				(size(sj, 1) == size(sj, 2)) || throw(DimensionMismatch("physical space mismatch"))			 						
 				s_p = space(sj, 1)
 				if isnothing(d)
 					d = s_p
 				else
 					(d == s_p) || throw(DimensionMismatch("physical space mismatch"))
 				end 				
			end
		end
	end
	isnothing(d) && throw(ArgumentError("pspace is missing"))
	for i in 1:m
		for j in 1:n
			sj = data[i, j]
			if isa(sj, AbstractMatrix)
				_is_id, scal = isid(sj)
				if _is_id
					sj = scal
				end
 			end
 			if isa(sj, AbstractMatrix)
 				sj = convert(M, sj)
 			else
 				isa(sj, Number) || throw(ArgumentError("elt should either be a tensor or a scalar"))
 				sj = convert(T, sj)
 			end
 			new_data[i, j] = sj
		end
	end
	return new_data, d
end
