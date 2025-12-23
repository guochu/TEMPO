function partialif_naive(lattice::AbstractADTLattice, rowind::ContourIndex, corr::AbstractCorrelationFunction, hyb::AdditiveHyb; 
						trunc::TruncationScheme=DefaultITruncation)
	ds = [lattice.d for i in 1:length(lattice)]
	b1 = branch(rowind)
	i = rowind.j
	i′ = (b1 == :τ) ? i+1 : i
	pos1 = index(lattice, i′, branch=b1)

	z = hyb.op
	(lattice.d == length(z)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	z2 = z .* z
	zz = reshape(kron(z, z), lattice.d, lattice.d)
	T = promote_type(scalartype(lattice), scalartype(hyb), scalartype(corr))
	tmp = vacuumstate(T, lattice)
	orth = Orthogonalize(SVD(), trunc)
	for b2 in branches(lattice)
		k2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
		for j in 1:k2
			j′ = (b1==:τ) ? j+1 : j
			pos2 = index(lattice, j′, branch=b2)
			coef = index(corr, i, j, b1=b1, b2=b2)

			if pos1 == pos2
				m = exp.(coef .* z2)
				t = ADTTerm((pos1, ), (m, ))
			else
				m = exp.(coef .* zz)
				t = ADTTerm((pos1, pos2), m)
			end
			apply!(t, tmp)
			canonicalize!(tmp, alg=orth)
		end
	end	
	return tmp
end

function partialif(lattice::AbstractADTLattice, rowind::ContourIndex, corr::AbstractCorrelationFunction, hyb::AdditiveHyb)
	ds = [lattice.d for i in 1:length(lattice)]
	b1 = branch(rowind)
	i = rowind.j
	i′ = (b1 == :τ) ? i+1 : i
	pos1 = index(lattice, i′, branch=b1)
	pos2s = Int[]
	coefs = scalartype(lattice)[]
	for b2 in branches(lattice)
		k2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
		for j in 1:k2
			j′ = (b1==:τ) ? j+1 : j
			pos2 = index(lattice, j′, branch=b2)
			coef = index(corr, i, j, b1=b1, b2=b2)
			push!(pos2s, pos2)
			push!(coefs, coef)			
		end
	end
	return partialif_densemps(ds, pos1, pos2s, hyb.op, coefs)
end

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
	return ADT(r)
end


