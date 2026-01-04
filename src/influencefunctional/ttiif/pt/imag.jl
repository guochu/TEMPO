# imaginary-time

"""
	influenceoperator(lattice::ImagGrassmannLattice1Order{<:A1Ā1B1B̄1}, corr2::ImagCorrelationFunction; algexpan)

The influenceoperator ΣᵢⱼΔᵢⱼāᵢaⱼ as an MPO, the bond dimension of MPO is 2n, where n is number of the prony expansion
"""
function influenceoperator(lattice::ImagPTLattice1Order, corr2::ImagCorrelationFunction, hyb::GeneralHybStyle; algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = transpose(corr2.data)
	op1, op2 = pairop(hyb)
	mpoj = pt_ti_mpotensor(corr, op1, op2, algexpan)
	mpotensors = _get_mpo3(mpoj)
	# println(mpotensors[2])
	return _fit_to_lattice(lattice, mpotensors) 
end

function influenceoperatorexponential(lattice::ImagPTLattice1Order, corr2::ImagCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::FirstOrderStepper; 
										algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = transpose(corr2.data)
	op1, op2 = pairop(hyb)
	mpoj = pt_ti_mpotensor(corr, op1, op2, algexpan)
	mpoj′ = timeevompo(mpoj, dt, alg)
	mpotensors = _get_mpo3(mpoj′)
	return _fit_to_lattice(lattice, mpotensors) 
end
function influenceoperatorexponential(lattice::ImagPTLattice1Order, corr2::ImagCorrelationFunction, dt::Real, hyb::GeneralHybStyle, alg::ComplexStepper; 
										algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = transpose(corr2.data)
	op1, op2 = pairop(hyb)
	mpoj = pt_ti_mpotensor(corr, op1, op2, algexpan)
	mpoja, mpojb = timeevompo(mpoj, dt, alg)
	mpo1, mpo2 = _get_mpo3(mpoja), _get_mpo3(mpojb)
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




function _get_mpo3(mpoj)
	# mpoj = ti_mpotensor(η, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	return tompotensors(h)
end

function pt_ti_mpotensor(corr::CorrelationMatrix, op1::AbstractMatrix, op2::AbstractMatrix, alg::ExponentialExpansionAlgorithm)
	# m1 = GenericDecayTerm(op1, op2, corr.ηₖⱼ[2:end])
	# m2 = GenericDecayTerm(op2, op1, corr.ηⱼₖ[2:end])
	m1 = GenericDecayTerm(op1, op2, corr.ηⱼₖ[2:end])
	m2 = GenericDecayTerm(op2, op1, corr.ηₖⱼ[2:end])


	m1s = exponential_expansion(m1, alg=alg)
	m2s = exponential_expansion(m2, alg=alg)

	# println("here---", corr.ηₖⱼ[1], " ", corr.ηⱼₖ[1])
	eta = corr.ηₖⱼ[1] + corr.ηⱼₖ[1]
	# h1 = corr.ηₖⱼ[1] * op1 * op2 + corr.ηⱼₖ[1] * op2 * op1
	# h1 = corr.ηₖⱼ[1] * op2 * op1 + corr.ηⱼₖ[1] * op1 * op2

	# h1 = (corr.ηₖⱼ[1] + corr.ηⱼₖ[1])  * op1 * op2 
	h1 = (eta/2) .* (op1 * op2 + op2 * op1)
	return SchurMPOTensor(h1, vcat(m1s, m2s))
	# return SchurMPOTensor(h1, [])
end
