abstract type AbstractFockOperator end

struct FockOperator{N, T <: Number} <: AbstractFockOperator
	positions::NTuple{N, Int}
	data::Vector{Array{T, 4}}
end


function FockOperator(positions::NTuple{N, Int}, data::NTuple{N, <:AbstractMatrix{T}}) where {N, T<:Number}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	p = TupleTools.sortperm(positions)
	positions = TupleTools.getindices(positions, p)
	data = TupleTools.getindices(data, p)
	return FockOperator(positions, [reshape(x, 1, size(x, 1), 1, size(x, 2)) for x in data])
end

TO.scalartype(::Type{FockOperator{N, T}}) where {N, T} = T

function apply!(x::FockTerm{N, T}, mps::FockMPO) where {N, T}
	data = x.data
	pos = x.positions
	pos_first = pos[1]
	pos_last = pos[end]

	ml = dropdims(data[1], dims=1)
	@tensor tmp[1,2,4,3,5] := mps[pos_first][1,2,3] * ml[4,5]
	mps[pos_first] = tie(n_fuse(tmp, 2), (1,1,2))
	if N == 1
		return mps
	end
	# println(size(mps[pos_first]))
	leftspace = space_r(data[1])
	for j in pos_first+1:pos_last-1
		posj = findfirst(y->y==j, pos)
		if isnothing(posj)
			d = size(mps[j], 2) 
			I2 = _eye(T, leftspace)
			mj = zeros(T, leftspace, d, leftspace)
			for i in 1:d
				mj[:, i, :] = I2
			end
		else
			mj = data[posj]
			leftspace = space_r(mj)
		end
		@tensor tmp[1,4,2,5,3,6] := mps[j][1,2,3] * mj[4,5,6]
		mps[j] = n_fuse(tie(tmp, (2,1,1,2)), 2)
	end
	mr = dropdims(data[end], dims=3)
	@tensor tmp[1,4,2,5,3] := mps[pos_last][1,2,3] * mr[4,5]
	mps[pos_last] = tie(n_fuse(tmp, 3), (2,1,1))
	for site in pos_first:pos_last-1
		@assert space_r(mps[site]) == space_l(mps[site+1])
	end
	return mps
end