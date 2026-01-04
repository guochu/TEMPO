# imaginary-time

"""
	influenceoperator(lattice::ImagGrassmannLattice1Order{<:A1Ā1B1B̄1}, corr2::ImagCorrelationFunction; algexpan)

The influenceoperator ΣᵢⱼΔᵢⱼāᵢaⱼ as an MPO, the bond dimension of MPO is 2n, where n is number of the prony expansion
"""
function influenceoperator(lattice::ImagPTLattice1Order, corr2::ImagCorrelationFunction, hyb::GeneralHybStyle; algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr2.data
	op1, op2 = pairop(hyb)
	mpoj = pt_ti_mpotensor(corr, op1, op2, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	mpotensors = tompotensors(h)
	# println(mpotensors[2])
	return _fit_to_lattice(lattice, mpotensors) 
end

function influenceoperatorexponential(lattice::ImagPTLattice1Order, corr2::ImagCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::FirstOrderStepper; 
										algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr2.data
	op1, op2 = pairop(hyb)
	mpoj = pt_ti_mpotensor(corr, op1, op2, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	h2 = timeevompo(h, dt, alg)
	mpotensors = tompotensors(h2)
	return _fit_to_lattice(lattice, mpotensors) 
end
function influenceoperatorexponential(lattice::ImagPTLattice1Order, corr2::ImagCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::ComplexStepper; 
										algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr2.data
	op1, op2 = pairop(hyb)
	mpoj = pt_ti_mpotensor(corr, op1, op2, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	h1, h2 = timeevompo(h, dt, alg)
	mpo1, mpo2 = tompotensors(h1), tompotensors(h2)
	return _fit_to_lattice(lattice, mpo1), _fit_to_lattice(lattice, mpo2) 
end


function differentialinfluencefunctional(lattice::ImagPTLattice1Order, corr::ImagCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::FirstOrderStepper, 
											algmult::DMRGAlgorithm; 
											algexpan::ExponentialExpansionAlgorithm=PronyExpansion()) 
	return influenceoperatorexponential(lattice, corr, dt, hyb, alg; algexpan=algexpan)
end
function differentialinfluencefunctional(lattice::ImagPTLattice1Order, corr::ImagCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::ComplexStepper, 
											algmult::DMRGAlgorithm; 
											algexpan::ExponentialExpansionAlgorithm=PronyExpansion()) 
	mpo1, mpo2 = influenceoperatorexponential(lattice, corr, dt, hyb, alg, algexpan=algexpan)
	return mult(mpo1, mpo2, algmult)
end

function _fit_to_lattice(lattice::ImagPTLattice1Order, mpotensors)
	L = length(lattice)
	data = similar(mpotensors, L)
	data[1] = mpotensors[1]
	data[end] = mpotensors[3]
	for j in 2:L-1
		data[j] = mpotensors[2]
	end
	return ProcessTensor(data)
end


# real-time
