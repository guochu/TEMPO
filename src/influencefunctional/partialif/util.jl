

# the algorithm in StrathearnLovett2018
function partialif_densemps(ds::Vector{Int}, row::Int, cols::Vector{Int}, op::Vector{<:Number}, coefs::Vector{<:Number})
	# println("row=", row, " cols ", cols)
	@assert length(cols) == length(coefs) > 1
	# @assert (row in cols)
	# @assert issorted(cols)
	p = sortperm(cols)
	cols = cols[p]
	coefs = coefs[p]
	T = promote_type(scalartype(op), scalartype(coefs)) 

	d = length(op)
	op2 = op .* op
	opop = reshape(kron(op, op), d, d)
	function onebody(m) 
		return exp.(m .* op2)
	end
	function twobody(m)
		return exp.(m .* opop)
	end
	L = length(cols)
	mpsdata = Vector{Array{T, 3}}(undef, L)
	pos = findfirst(x->row==x, cols)
	isnothing(pos) && throw(ArgumentError("$(row) is not a member of $cols"))
	# println("row = ", row)
	# println("cols = ", cols)
	if pos == 1
		m = onebody(coefs[1])
		tmp = zeros(T, d, d)
		for i in 1:d
			tmp[i, i] = m[i]
		end
		mpsdata[1] = reshape(tmp, 1,d,d)
		for j in 2:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, d,d,d)
			for i in 1:d
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp
		end
		m = twobody(coefs[L])
		mpsdata[L] = reshape(m,d,d,1)
	elseif pos == L
		m = twobody(coefs[1])
		mpsdata[1] = reshape(m, 1,d,d)
		for j in 2:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, d,d,d)
			for i in 1:d
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp			
		end
		m = onebody(coefs[L])
		tmp = zeros(T, d, d)
		for i in 1:d
			tmp[i, i] = m[i]
		end
		mpsdata[L] = reshape(tmp,d,d,1)
	else
		m = twobody(coefs[1])
		mpsdata[1] = reshape(m, 1,d,d)

		for j in 2:pos-1
			m = twobody(coefs[j])
			tmp = zeros(T, d,d,d)
			for i in 1:d
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp
		end

		m = onebody(coefs[pos])
		tmp = zeros(T,d,d,d)
		for i in 1:d
			tmp[i,i,i] = m[i]
		end
		mpsdata[pos] = tmp

		for j in pos+1:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, d,d,d)
			for i in 1:d
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp			
		end

		m = twobody(coefs[L])
		mpsdata[L] = reshape(m,d,d,1)
	end
	# return mpsdata
	return _fit_to_full(ds, cols, mpsdata)
end

function _fit_to_full(ds::Vector{Int}, pos, mpsdata)
	L = length(ds)
	r = similar(mpsdata, L)
	for j in 1:pos[1]-1
		r[j] = ones(1,ds[j],1)
	end
	for j in pos[end]+1:L
		r[j] = ones(1,ds[j],1)
	end
	leftspace = space_r(mpsdata[1])
	for j in pos[1]:pos[end]
		posj = findfirst(x->x==j, pos)
		if isnothing(posj)
			d = ds[j]
			tmp = zeros(leftspace,d,leftspace)
			for i in 1:d
				tmp[:,i,:] = _eye(leftspace)
			end
			r[j] = tmp
		else
			r[j] = mpsdata[posj]
			leftspace = space_r(mpsdata[posj])
		end
	end
	return FockMPS(r)
end


