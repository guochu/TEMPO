include("imag.jl")
include("real.jl")

function hybriddynamics!(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb, alg::TranslationInvariantIF)
	if alg.fast
		mps = hybriddynamics(lattice, corr, hyb, alg)
		return mult!(gmps, mps, alg.algmult)
	else
		(alg.verbosity > 1) && println("Serial scheme using 2^$(alg.k)-1 multiplications")
		return _hybriddynamics_slow!(gmps, lattice, corr, hyb, alg)
	end
end

function hybriddynamics(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb, alg::TranslationInvariantIF)
	if alg.fast
		(alg.verbosity > 1) && println("Tree bipartition scheme using $(alg.k) multiplications")
		return _hybriddynamics_fast(lattice, corr, hyb, alg)
	else
		(alg.verbosity > 1) && println("Serial scheme using 2^$(alg.k) multiplications")
		return _hybriddynamics_slow(lattice, corr, hyb, alg)
	end
end

function _hybriddynamics_fast(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb, alg::TranslationInvariantIF)
	algmult = alg.algmult
	if alg.verbosity > 1
		t = @elapsed mps = differentialinfluencefunctional(lattice, corr, 1/2^(alg.k), hyb, alg.algevo, algmult, algexpan=alg.algexpan)
		println("building the initial MPS-IF takes $t seconds, bond dimension is ", bond_dimension(mps))
	else
		mps = differentialinfluencefunctional(lattice, corr, 1/2^(alg.k), hyb, alg.algevo, algmult, algexpan=alg.algexpan)
	end
	
	for i in 1:alg.k
		if alg.verbosity > 1
			t = @elapsed mps = mult(mps, mps, algmult)
			println("the $i-th iteration takes $t seconds, bond dimension is ", bond_dimension(mps))
		else
			mps = mult(mps, mps, algmult)
		end
	end	
	return mps
end

function _hybriddynamics_slow(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb, alg::TranslationInvariantIF)
	algmult = alg.algmult
	if alg.verbosity > 1
		t = @elapsed mps0 = differentialinfluencefunctional(lattice, corr, 1/2^(alg.k), hyb, alg.algevo, algmult, algexpan=alg.algexpan)
		println("building the initial MPS-IF takes $t seconds, bond dimension is ", bond_dimension(mps))
	else
		mps0 = differentialinfluencefunctional(lattice, corr, 1/2^(alg.k), hyb, alg.algevo, algmult, algexpan=alg.algexpan)
	end
	mps = mps0

	for i in 1:2^(alg.k)-1
		if alg.verbosity > 1
			t = @elapsed mps = mult(mps, mps0, algmult)
			println("the $i-th iteration takes $t seconds, bond dimension is ", bond_dimension(mps))
		else
			mps = mult(mps, mps0, algmult)
		end		
	end
	return mps
end

function _hybriddynamics_slow!(gmps, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb, alg::TranslationInvariantIF)
	algmult = alg.algmult
	if alg.verbosity > 1
		t = @elapsed mps_all = influenceoperatorexponential(lattice, corr, 1/2^(alg.k), hyb, alg.algevo, algexpan=alg.algexpan)
		println("building the initial MPS-IFs takes $t seconds ")
	else
		mps_all = influenceoperatorexponential(lattice, corr, 1/2^(alg.k), hyb, alg.algevo, algexpan=alg.algexpan)
	end

	for i in 1:2^(alg.k)
		if alg.verbosity > 1
			t = @elapsed begin
				for mps in mps_all
					mult!(gmps, mps, algmult)
				end
			end
			println("the $i-th iteration takes $t seconds, bond dimension is ", bond_dimension(gmps))
		else
			for mps in mps_all
				mult!(gmps, mps, algmult)
			end
		end		
	end
	return gmps
end