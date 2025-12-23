abstract type AbstractFockTerm end

struct FockTerm{N, T <: Number} <: AbstractFockTerm
	data::NTuple{N, Array{T, 4}}
	positions::NTuple{N, Int}
end

function FockTerm(positions::NTuple{N, Int}, data::AbstractArray{T, M}) where {N, T, M}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	(2*N == M) || throw(DimensionMismatch("n positions mismatch with rank of data"))
	p = TupleTools.sortperm(positions)
	positions = TupleTools.getindices(positions, p)
	data = permute(data, [p; p .+ N])
	return FockTerm(positions, data)
end



FockTerm(p::Int, data::AbstractMatrix) = FockTerm((p,), (data,))
FockTerm(positions::AbstractVector{Int}, data::AbstractVector{<:AbstractMatrix}) = FockTerm((positions...,), (data...,))

TO.scalartype(::Type{FockTerm{N, T}}) where {N, T} = T


apply!(x::FockTerm{N, T}, mps::ProcessTensor; aheads::Union{NTuple{N,Bool}, Bool}=ntuple(j->true, N)) where {N, T} = apply!(x, mps, aheads)


function apply!(x::FockTerm{N, T}, mps::ProcessTensor, aheads::NTuple{N,Bool}) where {N, T}
	for (pos, m, a) in zip(x.positions, x.data, aheads)
		if a
			@tensor tmp[1,2,3,5] := mps[pos][1,2,3,4] * m[4,5]
		else
			@tensor tmp[3,1,4,5] := m[1,2] * mps[pos][3,2,4,5]
		end
		mps[pos] = tmp
	end
	return mps
end


apply!(x::FockTerm{N, T}, mps::ProcessTensor, ahead::Bool) where {N, T} = apply!(x, mps, ntuple(j->ahead, N))




