struct ContourOperator{T<:Number}
	indices::Vector{ContourIndex}
	ops::Vector{Matrix{T}}
end

function ContourOperator(idx::AbstractVector{ContourIndex}, data::AbstractVector{<:AbstractMatrix{T}}) where {T}
	(length(Set(idx)) == length(idx)) || throw(ArgumentError("multiple op with the same ContourIndex not allowed"))
	(length(idx) == length(data)) || throw(DimensionMismatch("op size mismatch with number of ContourIndexes"))
	return ContourOperator(convert(Vector{ContourIndex}, idx), convert(Vector{Matrix{T}}, data))
end
ContourOperator(p::ContourIndex, data::AbstractMatrix) = ContourOperator([p], [data])

TO.scalartype(::Type{ContourOperator{T}}) where {T} = T

apply!(x::ContourOperator, lat::AbstractPTLattice, mps::ProcessTensor; aheads::Union{AbstractVector{Bool}, Bool}=true) = apply!(x, lat, mps, aheads)

# function apply!(x::ContourOperator, lat::AbstractPTLattice, mps::ProcessTensor, aheads::AbstractVector{Bool})
# 	for (ind, m, a) in zip(x.indices, x.ops, aheads)
# 		pos = lat[ind]
# 		m2 = m
# 		if branch(ind) == :-
# 			m2 = transpose(m)
# 		end
# 		if a
# 			@tensor tmp[1,2,3,5] := mps[pos][1,2,3,4] * m2[4,5]
# 		else
# 			@tensor tmp[3,1,4,5] := m2[1,2] * mps[pos][3,2,4,5]
# 		end
# 		mps[pos] = tmp
# 	end
# 	return mps
# end

# apply!(x::ContourOperator, lat::AbstractPTLattice, mps::ProcessTensor, ahead::Bool) = apply!(x, lat, mps, [ahead for i in 1:length(x.indices)])


apply!(x::ContourOperator, lat::AbstractPTLattice, mps::ProcessTensor, aheads::Union{AbstractVector{Bool}, Bool}) = apply!(tofockprodterm(x, lat), mps, aheads)

function tofockprodterm(x::ContourOperator, lat::AbstractPTLattice) 
	L = length(x.indices)
	pos = Vector{Int}(undef, L)
	data = Vector{Matrix{scalartype(x)}}(undef, L)
	for i in 1:L
		ind = x.indices[i]
		pos[i] = lat[ind]
		data[i] = (branch(ind) == :-) ? transpose(x.ops[i]) : x.ops[i]
	end
	return ProdFockTerm(pos, data)
end