abstract type AbstractADTTerm end

struct ADTTerm{N, T <: Number} <: AbstractADTTerm
	data::Vector{Array{T, 3}}
	positions::NTuple{N, Int}
end

function ADTTerm(positions::NTuple{N, Int}, data::AbstractArray{T, N}) where {N, T<:Number}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	# iseven(N) || throw(ArgumentError("even number of variables expected"))
	p = TupleTools.sortperm(positions)
	positions = TupleTools.getindices(positions, p)
	data = permute(data, p)
	return ADTTerm(decompose_to_mps(data), positions)
end
ADTTerm(p::Int, data::AbstractVector{<:Number}) = ADTTerm((p,), data)
ADTTerm(p::Int, data::DenseMPSTensor) = ADTTerm((p,), [data])
function ADTTerm(positions::NTuple{N, Int}, data::AbstractVector{<:AbstractArray{T, 3}}) where {N, T}
	if issorted(positions)	
		for i in 1:length(positions)-1
			(space_r(data[i]) == space_l(data[i+1])) || throw(DimensionMismatch("MPO Tensor auxiliary space mismatch"))
		end
		(space_l(data[1]) == space_r(data[end])) || throw(DimensionMismatch("MPO Tensor auxiliary space mismatch"))
		return ADTTerm(data, positions)
	else
		throw(ArgumentError("not implemented for this case"))
	end
end

function ADTTerm(positions::NTuple{N, Int}, data::NTuple{N, <:AbstractVector{T}}) where {N, T<:Number}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	p = TupleTools.sortperm(positions)
	positions = TupleTools.getindices(positions, p)
	data = TupleTools.getindices(data, p)
	return ADTTerm([reshape(x, 1, length(x), 1) for x in data], positions)
end

TO.scalartype(::Type{ADTTerm{N, T}}) where {N, T} = T

function apply!(x::ADTTerm{N, T}, mps::ADT) where {N, T}
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

function decompose_to_mps(m::AbstractArray{T, N}) where {T<:Number, N}
	data = Vector{Array{T, 3}}(undef, N)
	if N == 1
		data[1] = reshape(m, 1, length(m), 1)
		return data
	end
	workspace = T[]
	q, r = tqr!(copy(m), (1,), ntuple(i->i+1, N-1), workspace)
	m = r
	data[1] = reshape(q, 1, size(q)...)
	L = N
	for i in 2:N-1
		q, r = tqr!(m, (1,2), ntuple(i->i+2, L-2), workspace)
		data[i] = q
		L -= 1
		m = r
	end
	data[N] = reshape(r, size(r)..., 1)
	return data
end