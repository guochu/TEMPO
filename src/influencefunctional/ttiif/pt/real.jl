# # real-time
# function influenceoperator(lattice::RealPTLattice1Order, corr2::RealCorrelationFunction, hyb::GeneralHybStyle; algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
# 	corr = corr2.data
# 	op1, op2 = pairop(hyb)
# 	mpoj1 = pt_ti_mpotensor(corr[:+, :+], op1, op2, algexpan)
# 	mpoj2 = pt_ti_mpotensor(corr[:+, :-], op1, op2, algexpan)
# 	mpoj3 = pt_ti_mpotensor(corr[:-, :+], op1, op2, algexpan)
# 	mpoj4 = pt_ti_mpotensor(corr[:-, :-], op1, op2, algexpan)
# 	h1, h2, h3, h4 = _get_mpo3(mpoj1), _get_mpo3(mpoj2), _get_mpo3(mpoj3), _get_mpo3(mpoj4)
# 	return _fit_to_lattice(lattice, mpotensors) 
# end

# function influenceoperatorexponential(lattice::RealPTLattice1Order, corr2::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::FirstOrderStepper; 
# 										algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
# 	corr = corr2.data
# 	op1, op2 = pairop(hyb)
# 	mpoj = pt_ti_mpotensor(corr, op1, op2, algexpan)
# 	h = MPOHamiltonian([mpoj, mpoj, mpoj])
# 	h2 = timeevompo(h, dt, alg)
# 	mpotensors = tompotensors(h2)
# 	return _fit_to_lattice(lattice, mpotensors) 
# end
# function influenceoperatorexponential(lattice::RealPTLattice1Order, corr2::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::ComplexStepper; 
# 										algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
# 	corr = corr2.data
# 	op1, op2 = pairop(hyb)
# 	mpoj = pt_ti_mpotensor(corr, op1, op2, algexpan)
# 	h = MPOHamiltonian([mpoj, mpoj, mpoj])
# 	h1, h2 = timeevompo(h, dt, alg)
# 	mpo1, mpo2 = tompotensors(h1), tompotensors(h2)
# 	return _fit_to_lattice(lattice, mpo1), _fit_to_lattice(lattice, mpo2) 
# end

# function _fit_to_lattice(lattice::RealPTLattice1Order, mpotensors, f1::Symbol, f2::Symbol)
# 	L = length(lattice)
# 	data = similar(mpotensors, L)
# 	T = scalartype(lattice)
# 	I2 = _eye(T, lattice.d)
# 	leftspace = 1

# 	j = lattice.N
# 	pos1, pos2 = index(lattice, j, branch=f1), index(lattice, j, branch=f2)
# 	if pos1 > pos2 # this sign has already been taken care of
# 		pos1, pos2 = pos2, pos1
# 	end
# 	if f1 = :+
		
# 	end



# 	data[1] = mpotensors[1]
# 	data[end] = mpotensors[3]
# 	for j in 2:L-1
# 		data[j] = mpotensors[2]
# 	end
# 	return ProcessTensor(data)
# end