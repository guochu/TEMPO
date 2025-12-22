abstract type AbstractFockTerm end

struct FockTerm{N, T <: Number} <: AbstractFockTerm
	positions::NTuple{N, Int}
	data::Vector{Array{T, 2}}
end

function FockTerm(positions::NTuple{N, Int}, data::AbstractVector{<:AbstractMatrix}) where {N}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	(N == length(data)) || throw(DimensionMismatch("op size mismatch with number of sites"))
	p = TupleTools.sortperm(positions)
	positions = TupleTools.getindices(positions, p)
	data = [data[pj] for pj in p]
	return FockTerm(positions, data)
end
FockTerm(p::Int, data::AbstractMatrix) = FockTerm((p,), [data])

TO.scalartype(::Type{FockTerm{N, T}}) where {N, T} = T

function apply!(x::FockTerm{N, T}, mps::ProcessTensor) where {N, T}
	for (pos, m) in zip(x.positions, x.data)
		@tensor tmp[3,1,4,5] := m[1,2] * mps[pos][3,2,4,5]
		mps[pos] = tmp
	end
	return mps
end