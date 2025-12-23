
function partialif_densempo(row::Int, cols::Vector{Int}, op::Matrix{<:Number}, coefs::Vector{<:Number})
	# println("row=", row, " cols ", cols)
	@assert length(cols) == length(coefs) > 1
	# @assert (row in cols)
	# @assert issorted(cols)
	p = sortperm(cols)
	cols = cols[p]
	coefs = coefs[p]
	T = promote_type(scalartype(op), scalartype(coefs)) 

	d = size(op, 1)
	d2 = d * d
	op2 = op .* op
	dim2 = CartesianIndices((d, d))
	opop = kron(op, op)
	function onebody(m) 
		return exp(m .* op2)
	end
	function twobody(m)
		return exp(m .* opop)
	end
	L = length(cols)
	mpsdata = Vector{Array{T, 4}}(undef, L)
	pos = findfirst(x->row==x, cols)
	isnothing(pos) && throw(ArgumentError("$(row) is not a member of $cols"))
	# println("row = ", row)
	# println("cols = ", cols)
	if pos == 1
		m = onebody(coefs[1])
		tmp = zeros(T, 1, d, d2, d)
		for i in 1:d2
			ind = dim2[i]
			tmp[1, ind[1], i, ind[2]] = 1
		end
		mpsdata[1] = _apply2(m, tmp)
		for j in 2:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, d2,d,d2,d)
			for i in 1:d2, j in 1:d2
				ind = dim2[j]
				tmp[i,ind[1],i, ind[2]] = m[i,j]
			end
			mpsdata[j] = tmp
		end
		m = twobody(coefs[L])
		tmp = zeros(T, d2,d,1,d)
		for i in 1:d2, j in 1:d2
			ind = dim2[j]
			tmp[i, ind[1], 1, ind[2]] = m[i, j]
		end
		mpsdata[L] = tmp
	elseif pos == L
		m = twobody(coefs[1])
		tmp = zeros(T, 1, d, d2, d)
		for i in 1:d2, j in 1:d2
			ind = dim2[i]
			tmp[1, ind[1], i, ind[2]] = m[i, j]
		end
		mpsdata[1] = tmp
		for j in 2:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, d2,d,d2,d)
			for i in 1:d2, j in 1:d2
				ind = dim2[j]
				tmp[i,ind[1],i, ind[2]] = m[i,j]
			end
			mpsdata[j] = tmp
		end
		tmp = zeros(T, d2,d,1,d)
		for i in 1:d2
			ind = dim2[i]
			tmp[i, ind[1], 1, ind[2]] = 1
		end
		m = onebody(coefs[L])
		mpsdata[L] = _apply2(m, tmp)
	else
		m = twobody(coefs[1])
		tmp = zeros(T, 1, d, d2, d)
		for i in 1:d2, j in 1:d2
			ind = dim2[i]
			tmp[1, ind[1], i, ind[2]] = m[i, j]
		end
		mpsdata[1] = tmp

		for j in 2:pos-1
			m = twobody(coefs[j])
			tmp = zeros(T, d2,d,d2,d)
			for i in 1:d2, j in 1:d2
				ind = dim2[j]
				tmp[i,ind[1],i, ind[2]] = m[i,j]
			end
			mpsdata[j] = tmp
		end

		m = onebody(coefs[pos])
		tmp = zeros(T,d2,d,d2,2)
		for i in 1:d
			ind = dim2[i]
			tmp[i,ind[1],i, ind[2]] = 1
		end
		mpsdata[pos] = _apply2(m, tmp)

		for j in pos+1:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, d2,d,d2,d)
			for i in 1:d2, j in 1:d2
				ind = dim2[j]
				tmp[i,ind[1],i, ind[2]] = m[i,j]
			end
			mpsdata[j] = tmp
		end

		m = twobody(coefs[L])
		tmp = zeros(T, d2,d,1,d)
		for i in 1:d2, j in 1:d2
			ind = dim2[j]
			tmp[i, ind[1], 1, ind[2]] = m[i, j]
		end
		mpsdata[L] = tmp
	end
	# return mpsdata
	return FockTerm(cols, mpsdata) 
end

function _apply2(m::AbstractMatrix, tmp::AbstractArray{<:Number, 4})
	@tensor tmp2[3,1,4,5] := m[1,2] * tmp[3,2,4,5]
end