abstract type AbstractFockTerm end

struct FockTerm{N, T <: Number} <: AbstractFockTerm
	positions::NTuple{N, Int}
	data::NTuple{N, Array{T, 2}}
end

function FockTerm(positions::NTuple{N, Int}, data::NTuple{N, <:AbstractMatrix}) where {N}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	p = TupleTools.sortperm(positions)
	positions = TupleTools.getindices(positions, p)
	data = TupleTools.getindices(data, p)
	return FockTerm(positions, data)
end
FockTerm(p::Int, data::AbstractMatrix) = FockTerm((p,), (data,))

TO.scalartype(::Type{FockTerm{N, T}}) where {N, T} = T

function apply!(x::FockTerm{N, T}, mps::ProcessTensor; ahead::Bool=true) where {N, T}
	if ahead
		for (pos, m) in zip(x.positions, x.data)
			# @tensor tmp[3,1,4,5] := m[1,2] * mps[pos][3,2,4,5]
			@tensor tmp[1,2,3,5] := mps[pos][1,2,3,4] * m[4,5]
			mps[pos] = tmp
		end
	else
		for (pos, m) in zip(x.positions, x.data)
			@tensor tmp[3,1,4,5] := m[1,2] * mps[pos][3,2,4,5]
			mps[pos] = tmp
		end		
	end
	return mps
end