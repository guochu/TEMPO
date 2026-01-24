abstract type AbstractFockTerm{T<:Number} end
TO.scalartype(::Type{<:AbstractFockTerm{T}}) where {T} = T
# num_terms(t::AbstractFockTerm) = length(t.positions)

struct FockTermS{N, T <: Number} <: AbstractFockTerm{T}
	data::NTuple{N, Array{T, 4}}
	positions::NTuple{N, Int}
end

function FockTermS(positions::NTuple{N, Int}, data::AbstractArray{T, M}) where {N, T, M}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	(2*N == M) || throw(DimensionMismatch("n positions mismatch with rank of data"))
	p = TupleTools.sortperm(positions)
	positions = TupleTools.getindices(positions, p)
	p2 = Vector{Int}(undef, M)
	for i in 1:N
		p2[2i-1] = p[i]
		p2[2i] = p[i] + N
	end
	data = permute(data, p2)
	r = decompose_to_mpo(data)
	return FockTermS(ntuple(j->r[j], N), positions)
end

# function FockTermS(positions::NTuple{N, Int}, data::NTuple{N, <:AbstractMatrix{T}}) where {N, T<:Number}
# 	for i in 1:N-1
# 		(positions[i] < positions[i+1]) || throw(ArgumentError("positions must be sorted"))
# 		(space_r(data[i]) == space_l(data[i+1])) || throw(DimensionMismatch("MPO Tensor auxiliary space mismatch"))
# 	end
# 	return FockTermS(data, positions)
# end

function FockTermS(positions::NTuple{N, Int}, data::NTuple{N, <:AbstractMatrix{T}}) where {N, T<:Number}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	p = TupleTools.sortperm(positions)
	positions = TupleTools.getindices(positions, p)
	data = TupleTools.getindices(data, p)
	return FockTermS(ntuple(j->_to4(data[j]), N), positions)
end


FockTermS(p::Int, data::AbstractMatrix) = FockTermS((p,), (data,))


struct FockTerm{T <: Number} <: AbstractFockTerm{T}
	data::Vector{Array{T, 4}}
	positions::Vector{Int}
end

function FockTerm(positions::AbstractVector{Int}, data::AbstractVector{<:AbstractArray{T, 4}}) where {T}
	(length(positions) == length(data)) || throw(DimensionMismatch("n positions mismatch with n data"))
	N = length(positions)
	for i in 1:N-1
		(positions[i] < positions[i+1]) || throw(ArgumentError("positions must be sorted"))
		(space_r(data[i]) == space_l(data[i+1])) || throw(DimensionMismatch("MPO Tensor auxiliary space mismatch"))
	end
	return FockTerm(convert(Vector{Array{T, 4}}, data), positions)
end

function apply!(x::AbstractFockTerm{T}, mps::ProcessTensor) where {T}
	data = x.data
	pos = x.positions
	N = length(data)
	pos_first = pos[1]
	pos_last = pos[end]

	ml = dropdims(data[1], dims=1)
	# println(size(ml), " ", size(mps[pos_first]))
	@tensor tmp[4,1,2,5,6] := ml[1,2,3] * mps[pos_first][4,3,5,6]
	mps[pos_first] = tie(tmp, (1,1,2,1))
	if N == 1
		return mps
	end
	leftspace = space_r(data[1])
	for j in pos_first+1:pos_last-1
		posj = findfirst(y->y==j, pos)
		if isnothing(posj)
			d = size(mps[j], 2) 
			Ia = _eye(T, leftspace)
			Id = _eye(T, d)
			@tensor mj[1,3,2,4] := Ia[1,2] * Id[3,4] 
		else
			mj = data[posj]
			leftspace = space_r(mj)
		end
		@tensor tmp[1,5,2,3,6,7] := mj[1,2,3,4] * mps[j][5,4,6,7] 
		mps[j] = tie(tmp, (2,1,2,1))
	end
	mr = dropdims(data[end], dims=3)
	@tensor tmp[1,4,2,5,6] := mr[1,2,3] * mps[pos_last][4,3,5,6] 
	mps[pos_last] = tie(tmp, (2,1,1,1))
	for site in pos_first:pos_last-1
		# println(size(mps[site]), " ", size( mps[site+1]))
		@assert space_r(mps[site]) == space_l(mps[site+1])
	end
	return mps
end


struct ProdFockTerm{T <: Number}
	data::Vector{Matrix{T}}
	positions::Vector{Int}
end

function ProdFockTerm(positions::AbstractVector{Int}, data::AbstractVector{<:AbstractMatrix{T}}) where {T}
	(length(positions) == length(data)) || throw(ArgumentError("number of positions mismatch with number of ops"))
	(length(positions) == length(Set(positions))) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	p = sortperm(positions)
	positions = positions[p]
	data = data[p]
	return ProdFockTerm(data, positions)
end

ProdFockTerm(pos::Int, data::AbstractMatrix{<:Number}) = ProdFockTerm([pos], [data])

TO.scalartype(::Type{ProdFockTerm{T}}) where {T} = T

apply!(x::ProdFockTerm, mps::ProcessTensor; aheads::Union{AbstractVector{Bool}, Bool}=true) = apply!(x, mps, aheads)
apply!(x::ProdFockTerm, mps::ProcessTensor, ahead::Bool) = apply!(x, mps, [ahead for i in 1:length(x.positions)])
function apply!(x::ProdFockTerm, mps::ProcessTensor, aheads::AbstractVector{Bool}) 
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

function decompose_to_mpo(m::AbstractArray{T, N}) where {T<:Number, N}
	L2 = div(N, 2)
	data = Vector{Array{T, 4}}(undef, L2)
	if L2 == 1
		data[1] = _to4(data[1])
		return data
	end
	workspace = T[]
	q, r = tqr!(copy(m), (1,2), ntuple(i->i+2, N-2), workspace)
	m = r
	d1, d2, d3 = size(q)
	data[1] = permute(reshape(q, 1, d1, d2, d3), (1,2,4,3))
	L = N
	for i in 2:L2-1
		q, r = tqr!(m, (1,2,3), ntuple(i->i+3, L-4), workspace)
		data[i] = permute(q, (1,2,4,3))
		L -= 2
		m = r
	end
	d1, d2, d3 = size(r)
	data[L2] = permute(reshape(r, d1, d2, d3, 1), (1,2,4,3))
	return data
end

function _to4(m::AbstractMatrix)
	d1, d2 = size(m)
	return reshape(m, (1,d1,1,d2))
end