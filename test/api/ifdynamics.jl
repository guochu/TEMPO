# Influence-functional dynamics on ADT / PT lattices: consistency of hybriddynamics /
# hybriddynamics! / naive variants, and cross-checks of XTRGIF and TDVPIF against the
# PartialIF reference dynamics.

@testset "hybriddynamics vs naive: ADT" begin
	δt = 0.1
	N = 2
	β = 2
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β)
	lattice = ADTLattice(N=N, δt=δt, contour=:real)
	corr = correlationfunction(bath, lattice)
	hyb = AdditiveHyb(diag(pauli_z()))
	trunc = truncdimcutoff(D=30, ϵ=1.0e-12)
	m1 = hybriddynamics(lattice, corr, hyb, trunc=trunc)
	m2 = hybriddynamics_naive(lattice, corr, hyb, trunc=trunc)
	@test distance(m1, m2) / norm(m2) < 1.0e-3
	m3 = hybriddynamics_naive!(vacuumstate(ComplexF64, lattice), lattice, corr, hyb, trunc=trunc)
	@test distance(m2, m3) < 1.0e-12
end

@testset "hybriddynamics vs naive: PT" begin
	δτ = 0.1
	N = 2
	β = N * δτ
	tol = 1.0e-4
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β)
	lattice = PTLattice(N=N, δτ=δτ, d=2, contour=:imag)
	corr = correlationfunction(bath, lattice)
	hyb = NonAdditiveHyb((m = randn(2, 2); m + m'))
	if1 = hybriddynamics(lattice, corr, hyb, PartialIF())
	if2 = hybriddynamics_naive(lattice, corr, hyb, PartialIF())
	@test distance(if1, if2) / norm(if1) < tol
end

# candidate algorithms to cross-check against the PartialIF reference: (algorithm, rtol);
# TDVPIF advances the influence operator stepwise with step δ
imag_algs(trunc, δτ) = [(XTRGIF(k=5, fast=true), 1.0e-2),
                        (XTRGIF(k=5, fast=false), 1.0e-2),
                        (TDVPIF(trunc=trunc, δ=δτ), 5.0e-2)]
real_algs(trunc, δt) = [(XTRGIF(k=5, algevo=WII(), algmult=SVDCompression(trunc)), 1.0e-2),
                        (XTRGIF(k=5, algmult=DMRG1(trunc=trunc, initguess=:svd)), 1.0e-2),
                        (XTRGIF(k=5, algmult=DMRG1(trunc=trunc), fast=false), 1.0e-2),
                        (TDVPIF(trunc=trunc, δ=δt), 5.0e-2)]

# cross-check the candidate algorithms against the PartialIF reference dynamics,
# plus one in-place flow from the vacuum with the `inplace` algorithm
function if_algorithm_checks(lattice, corr, hyb, base_alg, ref_naive::Bool, algs, inplace)
	mpsI1 = ref_naive ?
		hybriddynamics_naive(lattice, corr, hyb, base_alg) :
		hybriddynamics(lattice, corr, hyb, base_alg)
	for (alg, rtol) in algs
		mpsI2 = hybriddynamics(lattice, corr, hyb, alg)
		@test distance(mpsI1, mpsI2) / norm(mpsI1) < rtol
	end
	T0 = promote_type(scalartype(lattice), scalartype(corr), scalartype(hyb))
	mpsI2 = hybriddynamics!(vacuumstate(T0, lattice), lattice, corr, hyb, inplace)
	@test distance(mpsI1, mpsI2) / norm(mpsI1) < 1.0e-2
	return nothing
end

@testset "XTRGIF and TDVPIF vs PartialIF: imaginary-time" begin
	N = 6
	δτ = 0.1
	β = N * δτ
	d = 2
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=50, ϵ=1.0e-6, add_back=0)
	base_alg = PartialIF(trunc=trunc)
	algs = imag_algs(trunc, δτ)
	inplace = XTRGIF(k=5, fast=false)

	# ADT lattice, diagonal hybridization
	lattice = ADTLattice(N=N, δτ=δτ, d=d, contour=:imag)
	hyb = AdditiveHyb(normalize!(randn(Float64, d)))
	@testset "ADT" begin
		if_algorithm_checks(lattice, correlationfunction(bath, lattice), hyb, base_alg, false, algs, inplace)
	end

	# PT lattice, general hybridizations
	lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)
	corr = correlationfunction(bath, lattice)
	for T in (Float64, ComplexF64)
		@testset "PT ($T)" begin
			if_algorithm_checks(lattice, corr, NonAdditiveHyb(_rand_ham(T, d)), base_alg, true, algs, inplace)
		end
	end

	# in-place merge of the influence operator into pure impurity dynamics is
	# equivalent to multiplying the independently built IF into it
	for name in ("ADT", "PT")
		lattice = name == "ADT" ?
			ADTLattice(N=N, δτ=δτ, d=d, contour=:imag) :
			PTLattice(N=N, δτ=δτ, d=d, contour=:imag)
		corr = correlationfunction(bath, lattice)
		hyb = name == "ADT" ?
			AdditiveHyb(normalize!(randn(Float64, d))) :
			NonAdditiveHyb(_rand_ham(ComplexF64, d))
		@testset "TDVPIF merge: $name" begin
			model = ImpurityHamiltonian(_rand_ham(name == "ADT" ? Float64 : ComplexF64, d))
			mps4 = hybriddynamics!(sysdynamics(lattice, model, trunc=trunc), lattice, corr, hyb, TDVPIF(trunc=trunc, δ=δτ))
			g0 = sysdynamics(lattice, model, trunc=trunc)
			mult!(g0, hybriddynamics(lattice, corr, hyb, TDVPIF(trunc=trunc, δ=δτ)), SVDCompression(trunc))
			@test distance(mps4, g0) / norm(g0) < 5.0e-2
		end
	end
end

@testset "XTRGIF and TDVPIF vs PartialIF: real-time" begin
	N = 3
	δt = 0.1
	β = 1
	d = 2
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-6, add_back=0)
	base_alg = PartialIF(trunc=trunc)
	algs = real_algs(trunc, δt)
	inplace = XTRGIF(k=5, algmult=DMRG1(trunc=trunc), fast=false)

	# ADT lattice, diagonal hybridization
	lattice = ADTLattice(N=N, δt=δt, contour=:real, d=d)
	hyb = AdditiveHyb(normalize!(randn(Float64, d)))
	@testset "ADT" begin
		if_algorithm_checks(lattice, correlationfunction(bath, lattice), hyb, base_alg, false, algs, inplace)
	end

	# PT lattice, general hybridization
	lattice = PTLattice(N=N, δt=δt, contour=:real, d=d)
	hyb = NonAdditiveHyb(_rand_ham(ComplexF64, d))
	@testset "PT" begin
		if_algorithm_checks(lattice, correlationfunction(bath, lattice), hyb, base_alg, true, algs, inplace)
	end
end
