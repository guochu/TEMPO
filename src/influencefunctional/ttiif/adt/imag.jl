

"""
	influenceoperator(lattice, corr, hyb; algexpan=PronyExpansion())

Construct the (translation-invariant) influence functional as an MPO. The ADT version supports `ImagADTLattice1Order`/`RealADTLattice1Order` lattices with `AdditiveHyb` coupling; the PT version supports `ImagPTLattice1Order`/`RealPTLattice1Order` lattices with `GeneralHybStyle` (e.g. `NonAdditiveHyb`, `NonDiagonalHyb`) coupling.

# Arguments
- `lattice`: contour lattice (ADT or PT).
- `corr`: bath correlation function.
- `hyb`: system-bath coupling.
- `algexpan::ExponentialExpansionAlgorithm=PronyExpansion()`: exponential (Prony) expansion algorithm for the bath correlation function.

# Returns
The influence functional MPO: a single `ADT`/`ProcessTensor` in imaginary time; a tuple of 4 branch MPOs ((+,+), (+,−), (−,+), (−,−)) in real time.
"""
function influenceoperator(lattice::ImagADTLattice1Order, corr2::ImagCorrelationFunction, hyb::AdditiveHyb; algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr2.data
	op1, op2 = pairop(hyb)
	mpoj = adt_ti_mpotensor(corr, op1, op2, algexpan)
	mpstensors = _tompsj.(_get_mpo3(mpoj))
	return _fit_to_lattice(lattice, mpstensors) 
end

"""
	influenceoperatorexponential(lattice, corr, dt, hyb, alg; algexpan=PronyExpansion())

Construct the influence-functional exponential-operator MPO of a single time step (width `dt`) and time-evolve it with `alg` (`FirstOrderStepper` or `ComplexStepper`). `ComplexStepper` additionally returns the MPOs before and after evolution, so the returned tuple is twice as long as for `FirstOrderStepper` (1/2 in imaginary time and 4/8 in real time).

# Arguments
- `lattice`: contour lattice.
- `corr`: bath correlation function.
- `dt::Real`: width of a single time step.
- `hyb`: system-bath coupling.
- `alg`: time-evolution algorithm (`FirstOrderStepper` or `ComplexStepper`).
- `algexpan::ExponentialExpansionAlgorithm=PronyExpansion()`: exponential expansion algorithm.

# Returns
Tuple of time-evolved influence functional MPOs.
"""
function influenceoperatorexponential(lattice::ImagADTLattice1Order, corr2::ImagCorrelationFunction, dt::Real, hyb::AdditiveHyb, alg::FirstOrderStepper; 
										algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr2.data
	op1, op2 = pairop(hyb)
	mpoj = adt_ti_mpotensor(corr, op1, op2, algexpan)
	mpoj′ = timeevompo(mpoj, dt, alg)
	mpstensors = _tompsj.(_get_mpo3(mpoj′))
	return (_fit_to_lattice(lattice, mpstensors), )
end
function influenceoperatorexponential(lattice::ImagADTLattice1Order, corr2::ImagCorrelationFunction, dt::Real, hyb::AdditiveHyb, alg::ComplexStepper; 
										algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr2.data
	op1, op2 = pairop(hyb)
	mpoj = adt_ti_mpotensor(corr, op1, op2, algexpan)
	mpoja, mpojb = timeevompo(mpoj, dt, alg)
	mpo1, mpo2 = _tompsj.(_get_mpo3(mpoja)), _tompsj.(_get_mpo3(mpojb))
	return _fit_to_lattice(lattice, mpo1), _fit_to_lattice(lattice, mpo2) 
end


"""
	differentialinfluencefunctional(lattice, corr, dt, hyb, alg, algmult; algexpan=PronyExpansion())

Construct the differential influence functional, i.e. the full influence functional of a single time step, returned as an MPO/MPS. In real time, the branch MPOs obtained from `influenceoperatorexponential` are multiplied in order (with `algmult` controlling multiplication and compression); in imaginary time, the evolved operator is returned directly.

# Arguments
- `lattice`: contour lattice.
- `corr`: bath correlation function.
- `dt::Real`: width of a single time step.
- `hyb`: system-bath coupling.
- `alg`: time-evolution algorithm (`FirstOrderStepper` or `ComplexStepper`).
- `algmult::DMRGAlgorithm`: MPO multiplication (compression) algorithm.
- `algexpan::ExponentialExpansionAlgorithm=PronyExpansion()`: exponential expansion algorithm.

# Returns
The differential influence functional (`ADT` or `ProcessTensor`).
"""
function differentialinfluencefunctional(lattice::ImagADTLattice1Order, corr::ImagCorrelationFunction, dt::Real, hyb::AdditiveHyb, alg::FirstOrderStepper, 
											algmult::DMRGAlgorithm; 
											algexpan::ExponentialExpansionAlgorithm=PronyExpansion()) 
	mpo1, = influenceoperatorexponential(lattice, corr, dt, hyb, alg; algexpan=algexpan)
	return mpo1
end
function differentialinfluencefunctional(lattice::ImagADTLattice1Order, corr::ImagCorrelationFunction, dt::Real, hyb::AdditiveHyb, alg::ComplexStepper, 
											algmult::DMRGAlgorithm; 
											algexpan::ExponentialExpansionAlgorithm=PronyExpansion()) 
	mpo1, mpo2 = influenceoperatorexponential(lattice, corr, dt, hyb, alg, algexpan=algexpan)
	return mult(mpo1, mpo2, algmult)
end


function _fit_to_lattice(lattice::ImagADTLattice1Order, mpstensors)
	L = length(lattice)
	data = similar(mpstensors, L)
	j = lattice.k
	pos = index(lattice, j)
	@assert pos == 1
	data[pos] = mpstensors[1]
	for j in lattice.k-1:-1:3
		pos = index(lattice, j)
		data[pos] = mpstensors[2]
	end
	j = 2
	pos = index(lattice, j)
	data[pos] = mpstensors[3]
	j = 1
	pos = index(lattice, j)
	@assert pos == L
	data[pos] = ones(1, phydim(lattice), 1)
	return ADT(data)
end


function _tompsj(mpoj)
	a = ones(eltype(mpoj), size(mpoj, 4))
	@tensor tmp[1,2,3] := mpoj[1,2,3,4] * a[4]
	return tmp
end

function adt_ti_mpotensor(corr::CorrelationMatrix, op1::AbstractMatrix, op2::AbstractMatrix, alg::ExponentialExpansionAlgorithm)
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
	h1 = eta .* (op1 * op2 )
	return SchurMPOTensor(h1, vcat(m1s, m2s))
	# return SchurMPOTensor(h1, [])
end
