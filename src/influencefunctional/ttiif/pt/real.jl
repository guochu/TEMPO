# real-time
function influenceoperator(lattice::RealPTLattice1Order, corr::RealCorrelationFunction, hyb::GeneralHybStyle; algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr)
	op1, op2 = pairop(hyb)
	mpoj1 = pt_ti_mpotensor(η⁺⁺, op1, op2, algexpan)
	mpoj2 = pt_ti_mpotensor(η⁺⁻, fused_op(op1, :+), fused_op(op2, :-), algexpan)
	mpoj3 = pt_ti_mpotensor(η⁻⁺, fused_op(op1, :-), fused_op(op2, :+), algexpan)
	mpoj4 = pt_ti_mpotensor(η⁻⁻, transpose(op1), transpose(op2), algexpan)
	h1, h2, h3, h4 = _get_mpo3(mpoj1), _get_mpo3(mpoj2), _get_mpo3(mpoj3), _get_mpo3(mpoj4)
	mpo1 = _fit_to_lattice_diag(lattice, h1, :+, :+)
	mpo2 = _fit_to_lattice_offdiag(lattice, h2, :+, :-) 
	mpo3 = _fit_to_lattice_offdiag(lattice, h3, :-, :+) 
	mpo4 = _fit_to_lattice_diag(lattice, h4, :-, :-) 
	return mpo1, mpo2, mpo3, mpo4
end

function influenceoperatorexponential(lattice::RealPTLattice1Order, corr::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::FirstOrderStepper; 
										algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr)
	op1, op2 = pairop(hyb)
	mpoj1 = pt_ti_mpotensor(η⁺⁺, op1, op2, algexpan)
	mpoj2 = pt_ti_mpotensor(η⁺⁻, fused_op(op1, :+), fused_op(op2, :-), algexpan)
	mpoj3 = pt_ti_mpotensor(η⁻⁺, fused_op(op1, :-), fused_op(op2, :+), algexpan)
	mpoj4 = pt_ti_mpotensor(η⁻⁻, transpose(op1), transpose(op2), algexpan)
	mpoj1, mpoj2, mpoj3, mpoj4 = timeevompo(mpoj1, dt, alg), timeevompo(mpoj2, dt, alg), timeevompo(mpoj3, dt, alg), timeevompo(mpoj4, dt, alg)
	h1, h2, h3, h4 = _get_mpo3(mpoj1), _get_mpo3(mpoj2), _get_mpo3(mpoj3), _get_mpo3(mpoj4)
	mpo1 = _fit_to_lattice_diag(lattice, h1, :+, :+)
	mpo2 = _fit_to_lattice_offdiag(lattice, h2, :+, :-) 
	mpo3 = _fit_to_lattice_offdiag(lattice, h3, :-, :+) 
	mpo4 = _fit_to_lattice_diag(lattice, h4, :-, :-) 
	return mpo1, mpo2, mpo3, mpo4
end

function influenceoperatorexponential(lattice::RealPTLattice1Order, corr::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::ComplexStepper; 
										algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr)
	op1, op2 = pairop(hyb)
	mpoj1 = pt_ti_mpotensor(η⁺⁺, op1, op2, algexpan)
	mpoj2 = pt_ti_mpotensor(η⁺⁻, fused_op(op1, :+), fused_op(op2, :-), algexpan)
	mpoj3 = pt_ti_mpotensor(η⁻⁺, fused_op(op1, :-), fused_op(op2, :+), algexpan)
	mpoj4 = pt_ti_mpotensor(η⁻⁻, transpose(op1), transpose(op2), algexpan)
	mpoj1a, mpoj1b = timeevompo(mpoj1, dt, alg)
	mpoj2a, mpoj2b = timeevompo(mpoj2, dt, alg)
	mpoj3a, mpoj3b = timeevompo(mpoj3, dt, alg)
	mpoj4a, mpoj4b = timeevompo(mpoj4, dt, alg)
	h1a, h1b = _get_mpo3(mpoj1a), _get_mpo3(mpoj1b)
	h2a, h2b = _get_mpo3(mpoj2a), _get_mpo3(mpoj2b)
	h3a, h3b = _get_mpo3(mpoj3a), _get_mpo3(mpoj3b)
	h4a, h4b = _get_mpo3(mpoj4a), _get_mpo3(mpoj4b)
	mpo1a, mpo1b = _fit_to_lattice_diag(lattice, h1a, :+, :+), _fit_to_lattice_diag(lattice, h1b, :+, :+)
	mpo2a, mpo2b = _fit_to_lattice_offdiag(lattice, h2a, :+, :-), _fit_to_lattice_offdiag(lattice, h2b, :+, :-) 
	mpo3a, mpo3b = _fit_to_lattice_offdiag(lattice, h3a, :-, :+), _fit_to_lattice_offdiag(lattice, h3b, :-, :+) 
	mpo4a, mpo4b = _fit_to_lattice_diag(lattice, h4a, :-, :-), _fit_to_lattice_diag(lattice, h4b, :-, :-) 
	return (mpo1a, mpo1b), (mpo2a, mpo2b), (mpo3a, mpo3b), (mpo4a, mpo4b)
end


function differentialinfluencefunctional(lattice::RealPTLattice1Order, corr::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::FirstOrderStepper, 
											algmult::DMRGAlgorithm; algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	h1, h2, h3, h4 = influenceoperatorexponential(lattice, corr, dt, hyb, alg, algexpan=algexpan)
	# trunc = algmult.trunc
	# canonicalize!(h1, alg=Orthogonalize(trunc=trunc, normalize=false))
	# canonicalize!(h2, alg=Orthogonalize(trunc=trunc, normalize=false))
	# canonicalize!(h3, alg=Orthogonalize(trunc=trunc, normalize=false))
	# canonicalize!(h4, alg=Orthogonalize(trunc=trunc, normalize=false))
	mps = mult(h2, h1, algmult)
	mps = mult(h3, mps, algmult)
	mps = mult(h4, mps, algmult)
	return mps
end
function differentialinfluencefunctional(lattice::RealPTLattice1Order, corr::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::ComplexStepper, 
											algmult::DMRGAlgorithm; algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	(h1a, h1b), (h2a, h2b), (h3a, h3b), (h4a, h4b) = influenceoperatorexponential(lattice, corr, dt, hyb, alg, algexpan=algexpan)
	# trunc = algmult.trunc
	# canonicalize!(h1a, alg=Orthogonalize(trunc=trunc, normalize=false))
	# canonicalize!(h1b, alg=Orthogonalize(trunc=trunc, normalize=false))
	# canonicalize!(h2a, alg=Orthogonalize(trunc=trunc, normalize=false))
	# canonicalize!(h2b, alg=Orthogonalize(trunc=trunc, normalize=false))
	# canonicalize!(h3a, alg=Orthogonalize(trunc=trunc, normalize=false))
	# canonicalize!(h3b, alg=Orthogonalize(trunc=trunc, normalize=false))
	# canonicalize!(h4a, alg=Orthogonalize(trunc=trunc, normalize=false))
	# canonicalize!(h4b, alg=Orthogonalize(trunc=trunc, normalize=false))

	mps = mult(h1b, h1a, algmult)

	mps = mult(h2a, mps, algmult)
	mps = mult(h2b, mps, algmult)

	mps = mult(h3a, mps, algmult)
	mps = mult(h3b, mps, algmult)

	mps = mult(h4a, mps, algmult)
	mps = mult(h4b, mps, algmult)

	return mps
end

function _get_signed_corr(lattice::RealPTLattice1Order, corr::RealCorrelationFunction)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = branch(corr, :+, :+), branch(corr, :+, :-), branch(corr, :-, :+), branch(corr, :-, :-)
	# if index(lattice, 1, branch=:+) > index(lattice, 1, branch=:+)
	# 	η⁺⁺ = transpose(η⁺⁺)
	# end
	# if index(lattice, 1, branch=:-) > index(lattice, 1, branch=:+)
	# 	η⁺⁻ = transpose(η⁺⁻)
	# end
	# if index(lattice, 1, branch=:+) > index(lattice, 1, branch=:-)
	# 	η⁻⁺ = transpose(η⁻⁺)
	# end
	# if index(lattice, 1, branch=:-) > index(lattice, 1, branch=:-)
	# 	η⁻⁻ = transpose(η⁻⁻)
	# end
	return η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻
end

# _contour_op(op1::AbstractMatrix, f1::Symbol) = (f1 == :-) ? transpose(op1) : op1

function fused_op(op1::AbstractMatrix, f1::Symbol)
	I2 = one(op1)
	d = size(I2, 1)
	d2 = length(I2)
	f = reshape(_eye(scalartype(I2), d2, d2), d2, d, d)
	if f1 == :+
		x = op1
		y = I2
	else
		x = I2
		y = transpose(op1)
	end
	@tensor a[5,6] := x[1,2] * y[3,4] * f[5,1,3] * f[6,2,4]
	return a
end

function _fit_to_lattice_diag(lattice::RealPTLattice1Order, mpotensors, f1::Symbol, f2::Symbol, trunc::TruncationScheme=DefaultMPOTruncation)
	@assert length(mpotensors) == 3
	L = length(lattice)
	data2 = similar(mpotensors, L)
	d = phydim(lattice)
	T = scalartype(lattice)
	I2 = _eye(T, d)

	leftspace = 1
	j = lattice.N
	pos1, pos2 = index(lattice, j, branch=:+), index(lattice, j, branch=:-)
	@assert pos1 + 1 == pos2 == 2
	if f1 == :+
		data2[pos1] = mpotensors[1]
		leftspace = space_r(data2[pos1])
		vd2 = _eye(T, leftspace)
		data2[pos2] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
	else
		vd2 = _eye(T, leftspace)
		data2[pos1] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
		data2[pos2] = mpotensors[1]
	end
	leftspace = space_r(data2[pos2])
	for j in lattice.N-1:-1:2
		pos1, pos2 = index(lattice, j, branch=:+), index(lattice, j, branch=:-)
		@assert pos1 + 1 == pos2
		if f1 == :+
			data2[pos1] = mpotensors[2]
			leftspace = space_r(data2[pos1])
			vd2 = _eye(T, leftspace)
			data2[pos2] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
		else
			vd2 = _eye(T, leftspace)
			data2[pos1] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
			data2[pos2] = mpotensors[2]
		end
		leftspace = space_r(data2[pos2])
	end
	j = 1
	pos1, pos2 = index(lattice, j, branch=:+), index(lattice, j, branch=:-)
	@assert pos1 + 1 == pos2 == L
	if f1 == :+
		data2[pos1] = mpotensors[3]
		leftspace = space_r(data2[pos1])
		vd2 = _eye(T, leftspace)
		data2[pos2] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
	else
		vd2 = _eye(T, leftspace)
		data2[pos1] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
		data2[pos2] = mpotensors[3]
	end

	# for i in 1:length(data2)-1
	# 	println("i = ", i, " ", size(data2[i]), " ", size(data2[i+1]))
	# end

	return ProcessTensor(data2)
end

function _fit_to_lattice_offdiag(lattice::RealPTLattice1Order, mpotensors, f1::Symbol, f2::Symbol, trunc::TruncationScheme=DefaultMPOTruncation)
	@assert length(mpotensors) == 3
	L = length(lattice)
	data2 = similar(mpotensors, L)
	d = phydim(lattice)
	T = scalartype(lattice)
	I2 = _eye(T, d)

	u_left, v_left = split_mpotensor(mpotensors[1], trunc)
	u_middle, v_middle = split_mpotensor(mpotensors[2], trunc)
	u_right, v_right = split_mpotensor(mpotensors[3], trunc)
	# println(size(u_left), " ", size(v_left), " ", size(u_middle), " ", size(v_middle), " ", size(u_right), " ", size(v_right))

	leftspace = 1

	j = lattice.N
	pos1, pos2 = index(lattice, j, branch=f1), index(lattice, j, branch=f2)
	if pos1 > pos2 # this sign has already been taken care of
		# println("here----xxxxx", f1, " ", f2)
		pos1, pos2 = pos2, pos1
	end
	posa, posb = band_boundary(lattice, j)
	# println(posa, " ", posb, " ", pos1, " ", pos2)

	for i in posa:posb
		if i < pos1
			vd2 = _eye(T, leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
		elseif i == pos1
			data2[i] = u_left
		elseif i < pos2
			vd2 = _eye(T, leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
		elseif i == pos2
			data2[i] = v_left
		else
			vd2 = _eye(T, leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
		end
		leftspace = space_r(data2[i])
	end
	for j in lattice.N-1:-1:2
		posa, posb = band_boundary(lattice, j)
		pos1, pos2 = index(lattice, j, branch=f1), index(lattice, j, branch=f2)
		if pos1 > pos2 # this sign has already been taken care of
			pos1, pos2 = pos2, pos1
		end
		for i in posa:posb
			if i < pos1
				vd2 = _eye(T, leftspace)
				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
			elseif i == pos1
				data2[i] = u_middle
			elseif i < pos2
				vd2 = _eye(T, leftspace)
				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
			elseif i == pos2
				data2[i] = v_middle
			else
				vd2 = _eye(T, leftspace)
				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
			end
			leftspace = space_r(data2[i])
		end
	end
	j = 1
	posa, posb = band_boundary(lattice, j)
	pos1, pos2 = index(lattice, j, branch=f1), index(lattice, j, branch=f2)
	if pos1 > pos2 # this sign has already been taken care of
		pos1, pos2 = pos2, pos1
	end
	for i in posa:posb
		if i < pos1
			vd2 = _eye(T, leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
		elseif i == pos1
			data2[i] = u_right
		elseif i < pos2
			vd2 = _eye(T, leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
		elseif i == pos2
			data2[i] = v_right
		else
			vd2 = _eye(T, leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
		end
		leftspace = space_r(data2[i])
	end

	# for i in 1:length(data2)-1
	# 	println("i = ", i, " ", size(data2[i]), " ", size(data2[i+1]))
	# end

	return ProcessTensor(data2)
end

band_boundary(lattice::RealPTLattice1Order, j::Int) = index(lattice, j, branch=:+), index(lattice, j, branch=:-)


function split_mpotensor(mpoj::DenseMPOTensor, trunc)
	d2 = size(mpoj, 2)
	# println("here---", size(mpoj))
	d = round(Int, sqrt(d2))
	@assert d*d == d2
	f = reshape(_eye(scalartype(mpoj), d2, d2), d2, d, d)
	@tensor mpoj6[1,5,7, 6,3,8] := mpoj[1,2,3,4] * f[2,5,6] * f[4,7,8]
	u, s, v = tsvd!(mpoj6, (1,2,3), (4,5,6), trunc=trunc)
	ss = Matrix(Diagonal(sqrt.(s)))
	@tensor uu[1,2,5,3] := u[1,2,3,4] * ss[4,5]
	@tensor vv[1,3,4,5] := ss[1,2] * v[2,3,4,5]
	return uu, vv
end
