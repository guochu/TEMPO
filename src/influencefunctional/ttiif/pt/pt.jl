include("imag.jl")
include("real.jl")


function hybriddynamics!(gmps::ProcessTensor, lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, hyb::GeneralHybStyle, alg::TranslationInvariantIF)
	mps = hybriddynamics(lattice, corr, hyb, alg)
	return mult!(gmps, mps, alg.algmult)
end

function hybriddynamics(lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, hyb::GeneralHybStyle, alg::TranslationInvariantIF)
	if alg.fast
		(alg.verbosity > 1) && println("Tree bipartition scheme using $(alg.k) multiplications")
		return _hybriddynamics_fast(lattice, corr, hyb, alg)
	else
		(alg.verbosity > 1) && println("Serial scheme using 2^$(alg.k)-1 multiplications")
		return _hybriddynamics_slow(lattice, corr, hyb, alg)
	end
end

function _hybriddynamics_fast(lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, hyb::GeneralHybStyle, alg::TranslationInvariantIF)
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

function _hybriddynamics_slow(lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, hyb::GeneralHybStyle, alg::TranslationInvariantIF)
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