function partialif_naive(lattice::AbstractPTLattice, rowind::ContourIndex, corr::AbstractCorrelationFunction, hyb::NonAdditiveHyb; 
							trunc::TruncationScheme=DefaultITruncation)
	z = hyb.op
	(lattice.d == size(z, 1) == size(z, 2)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	# d = lattice.d
	# z2 = z * z
	# zz = kron(z, z)

	b1 = branch(rowind)
	i = rowind.j
	# pos1 = index(lattice, i, branch=b1)
		
	T = promote_type(scalartype(lattice), scalartype(hyb), scalartype(corr))
	tmp = vacuumstate(T, lattice)
	orth = Orthogonalize(SVD(), trunc)
	for b2 in branches(lattice)
		k2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
		for j in 1:k2
			coef = index(corr, i, j, b1=b1, b2=b2)
			ind2 = ContourIndex(j, b2) 
			# pos2 = lattice[ind2]
			# if pos1 == pos2
			# 	m = exp(coef .* z2)
			# 	# t = ContourOperator(ind1, m)
			# 	t = FockTermS(pos1, m)
			# else
			# 	m = exp(coef .* zz)
			# 	t = FockTermS((pos1, pos2), reshape(m, (d,d,d,d)))
			# end
			t = _get_contour_op(lattice, rowind, ind2, z, coef)
			apply!(t, tmp)
			canonicalize!(tmp, alg=orth)
		end
	end
	return tmp
end

function _get_contour_op(lattice, ind1::ContourIndex, ind2::ContourIndex, z::AbstractMatrix, coef)
	d = lattice.d
	z1 = z
	z2 = z
	if branch(ind1) == :-
		z1 = transpose(z)
	end
	if branch(ind2) == :-
		z2 = transpose(z)
	end
	pos1, pos2 = lattice[ind1], lattice[ind2]
	if pos1 == pos2
		m = exp(coef .* z1 * z2) 
		t = FockTermS(pos1, m)
	else
		zz = kron(z2, z1)
		m = exp(coef .* zz)
		t = FockTermS((pos1, pos2), reshape(m, (d,d,d,d)))
	end
	return t
end

# function partialif(lattice::AbstractPTLattice, rowind::ContourIndex, corr::AbstractCorrelationFunction, hyb::NonAdditiveHyb)
# 	k = lattice.N
# 	(lattice.d == size(hyb.op, 1)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))

# 	b1 = branch(rowind)
# 	i = rowind.j
# 	pos1 = index(lattice, i, branch=b1)

# 	pos2s = Int[]
# 	coefs = scalartype(lattice)[]

# 	for b2 in branches(lattice)
# 		k2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
# 		for j in 1:k2
# 			pos2 = index(lattice, j, branch=b2)
# 			coef = index(corr, i, j, b1=b1, b2=b2)
# 			push!(pos2s, pos2)
# 			push!(coefs, coef)			
# 		end
# 	end
# 	pos2s, mpsdata = partialif_densempo(pos1, pos2s, hyb.op, coefs)
# 	return _fit_to_full(length(lattice), lattice.d, pos2s, mpsdata)
# end


# function partialif_densempo(row::Int, cols::Vector{Int}, op::Matrix{<:Number}, coefs::Vector{<:Number})
# 	# println("row=", row, " cols ", cols)
# 	@assert length(cols) == length(coefs) > 1
# 	# @assert (row in cols)
# 	# @assert issorted(cols)
# 	p = sortperm(cols)
# 	cols = cols[p]
# 	coefs = coefs[p]
# 	T = promote_type(scalartype(op), scalartype(coefs)) 

# 	d = size(op, 1)
# 	d2 = d * d
# 	op2 = op * op
# 	dim2 = CartesianIndices((d, d))
# 	opop = kron(op, op)
# 	function onebody(m) 
# 		return exp(m .* op2)
# 	end
# 	function twobody(m)
# 		# return exp(m .* opop)
# 		m2 = exp(m .* opop)
# 		# return m2
# 		return reshape(permute(reshape(m2, d, d, d, d), (1,3,2,4)), d2, d2)
# 	end
# 	L = length(cols)
# 	mpsdata = Vector{Array{T, 4}}(undef, L)
# 	pos = findfirst(x->row==x, cols)
# 	isnothing(pos) && throw(ArgumentError("$(row) is not a member of $cols"))
# 	# println("row = ", row, " pos = ", pos, " L = ", L)
# 	# println("cols = ", cols, " coefs = ", coefs)
# 	if pos == 1
# 		println("1------------------------")
# 		m = onebody(coefs[1])
# 		tmp = zeros(T, 1, d, d2, d)
# 		for i in 1:d2
# 			ind = dim2[i]
# 			tmp[1, ind[1], i, ind[2]] = 1
# 		end
# 		mpsdata[1] = _apply2(m, tmp)
# 		for j in 2:L-1
# 			m = twobody(coefs[j])
# 			tmp = zeros(T, d2,d,d2,d)
# 			for i in 1:d2, k in 1:d2
# 				ind = dim2[k]
# 				tmp[i,ind[1],i, ind[2]] = m[i,k]
# 			end
# 			mpsdata[j] = tmp
# 		end
# 		m = twobody(coefs[L])
# 		tmp = zeros(T, d2,d,1,d)
# 		for i in 1:d2, j in 1:d2
# 			ind = dim2[j]
# 			tmp[i, ind[1], 1, ind[2]] = m[i, j]
# 		end
# 		mpsdata[L] = tmp
# 	elseif pos == L
# 		println("2------------------------")
# 		m = twobody(coefs[1])
# 		tmp = zeros(T, 1, d, d2, d)
# 		for i in 1:d2, j in 1:d2
# 			ind = dim2[j]
# 			tmp[1, ind[1], i, ind[2]] = m[i, j]
# 		end
# 		mpsdata[1] = tmp
# 		for k in 2:L-1
# 			m = twobody(coefs[k])
# 			tmp = zeros(T, d2,d,d2,d)
# 			for i in 1:d2, j in 1:d2
# 				ind = dim2[j]
# 				tmp[i,ind[1],i, ind[2]] = m[i, j]
# 			end
# 			mpsdata[k] = tmp
# 		end
# 		tmp = zeros(T, d2,d,1,d)
# 		for i in 1:d2
# 			ind = dim2[i]
# 			tmp[i, ind[1], 1, ind[2]] = 1
# 		end
# 		m = onebody(coefs[L])
# 		mpsdata[L] = _apply2(m, tmp)
# 	else
# 		println("3------------------------")
# 		m = twobody(coefs[1])
# 		tmp = zeros(T, 1, d, d2, d)
# 		for i in 1:d2, j in 1:d2
# 			ind = dim2[i]
# 			tmp[1, ind[1], i, ind[2]] = m[i, j]
# 		end
# 		mpsdata[1] = tmp

# 		for j in 2:pos-1
# 			m = twobody(coefs[j])
# 			tmp = zeros(T, d2,d,d2,d)
# 			for i in 1:d2, k in 1:d2
# 				ind = dim2[k]
# 				tmp[i,ind[1],i, ind[2]] = m[i,k]
# 			end
# 			mpsdata[j] = tmp
# 		end

# 		m = onebody(coefs[pos])
# 		tmp = zeros(T,d2,d,d2,2)
# 		for i in 1:d
# 			ind = dim2[i]
# 			tmp[i,ind[1],i, ind[2]] = 1
# 		end
# 		mpsdata[pos] = _apply2(m, tmp)

# 		for j in pos+1:L-1
# 			m = twobody(coefs[j])
# 			tmp = zeros(T, d2,d,d2,d)
# 			for i in 1:d2, k in 1:d2
# 				ind = dim2[k]
# 				tmp[i,ind[1],i, ind[2]] = m[i,k]
# 			end
# 			mpsdata[j] = tmp
# 		end

# 		m = twobody(coefs[L])
# 		tmp = zeros(T, d2,d,1,d)
# 		for i in 1:d2, j in 1:d2
# 			ind = dim2[j]
# 			tmp[i, ind[1], 1, ind[2]] = m[i, j]
# 		end
# 		mpsdata[L] = tmp
# 	end
# 	# return mpsdata
# 	return cols, mpsdata
# end

# function _fit_to_full(L::Int, d::Int, pos, mpsdata)
# 	r = similar(mpsdata, L)
# 	I2 = _eye(d)
# 	for j in 1:pos[1]-1
# 		r[j] = reshape(I2, 1, d, 1, d)
# 	end
# 	for j in pos[end]+1:L
# 		r[j] = reshape(I2, 1, d, 1, d)
# 	end
# 	leftspace = space_r(mpsdata[1])
# 	for j in pos[1]:pos[end]
		
# 		posj = findfirst(x->x==j, pos)
# 		if isnothing(posj)
# 			Ia = _eye(leftspace, leftspace)
# 			@tensor mj[1,3,2,4] := Ia[1,2] * I2[3,4] 
# 		else
# 			mj = mpsdata[posj]
# 			leftspace = space_r(mj)
# 		end
# 		r[j] = mj
# 	end
# 	return ProcessTensor(r)
# end

# function _apply2(m::AbstractMatrix, tmp::AbstractArray{<:Number, 4})
# 	@tensor tmp2[3,1,4,5] := m[1,2] * tmp[3,2,4,5]
# end