println("------------------------------------")
println("|              TDVPIF              |")
println("------------------------------------")

@testset "TDVPIF: ADT imaginary-time" begin

	N = 6
	δτ = 0.1
	β = N * δτ
	rtol = 5.0e-2

	trunc = truncdimcutoff(D=50, ϵ=1.0e-10, add_back=0)

	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	corr = Δτ(bath, N=N, δτ=δτ)

	d = 2
	lattice = ADTLattice(N=N, δτ=δτ, d=d, contour=:imag)

	op = randn(Float64, d)
	op ./= norm(op)
	hyb = AdditiveHyb(op)

	mps1 = hybriddynamics(lattice, corr, hyb, PartialIF(trunc=trunc))

	for δ in (0.1, 0.05)
		mps2 = hybriddynamics(lattice, corr, hyb, TDVPIF(trunc=trunc, δ=δ))
		@test distance(mps1, mps2) / norm(mps1) < rtol
	end

	mps2 = hybriddynamics!(vacuumstate(Float64, lattice), lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05))
	@test distance(mps1, mps2) / norm(mps1) < rtol

	# in-place flow on a non-trivial initial MPO: the influence operator is
	# merged into the pure impurity dynamics (from sysdynamics) in one flow,
	# equivalent to multiplying the independently built IF into it
	model = ImpurityHamiltonian(_rand_ham(Float64, d))
	mps4 = hybriddynamics!(sysdynamics(lattice, model, trunc=trunc), lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05))
	g0 = sysdynamics(lattice, model, trunc=trunc)
	mult!(g0, hybriddynamics(lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05)), SVDCompression(trunc))
	@test distance(mps4, g0) / norm(g0) < rtol
end


@testset "TDVPIF: PT imaginary-time" begin

	N = 6
	δτ = 0.1
	β = N * δτ
	rtol = 5.0e-2

	trunc = truncdimcutoff(D=50, ϵ=1.0e-10, add_back=0)

	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	corr = Δτ(bath, N=N, δτ=δτ)

	d = 2
	lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)

	for T in (Float64, ComplexF64)
		hyb = NonAdditiveHyb(_rand_ham(T, d))

		mps1 = hybriddynamics_naive(lattice, corr, hyb, PartialIF(trunc=trunc))

		for δ in (0.1, 0.05)
			mps2 = hybriddynamics(lattice, corr, hyb, TDVPIF(trunc=trunc, δ=δ))
			@test distance(mps1, mps2) / norm(mps1) < rtol
		end

		mps2 = hybriddynamics!(vacuumstate(promote_type(scalartype(lattice), scalartype(corr), T), lattice), lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05))
		@test distance(mps1, mps2) / norm(mps1) < rtol

		# in-place flow on a non-trivial initial MPO (pure impurity dynamics)
		model = ImpurityHamiltonian(_rand_ham(T, d))
		mps4 = hybriddynamics!(sysdynamics(lattice, model, trunc=trunc), lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05))
		g0 = sysdynamics(lattice, model, trunc=trunc)
		mult!(g0, hybriddynamics(lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05)), SVDCompression(trunc))
		@test distance(mps4, g0) / norm(g0) < rtol
	end
end


@testset "TDVPIF: ADT real-time" begin

	N = 3
	δt = 0.1
	β = 1

	rtol = 5.0e-2
	trunc = truncdimcutoff(D=100, ϵ=1.0e-6, add_back=0)

	spec = Leggett(d=1, ωc=1)

	bath = bosonicbath(spec, β=β, μ=0)
	corr = Δt(bath, N=N, t=N*δt)

	d = 2
	lattice = ADTLattice(N=N, δt=δt, contour=:real, d=d)

	op = randn(Float64, d)
	op ./= norm(op)
	hyb = AdditiveHyb(op)

	mps1 = hybriddynamics(lattice, corr, hyb, PartialIF(trunc=trunc))

	for δ in (0.1, 0.05)
		mps2 = hybriddynamics(lattice, corr, hyb, TDVPIF(trunc=trunc, δ=δ))
		@test distance(mps1, mps2) / norm(mps1) < rtol
	end

	mps2 = hybriddynamics!(vacuumstate(ComplexF64, lattice), lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05))
	@test distance(mps1, mps2) / norm(mps1) < rtol

	# in-place flow on a non-trivial initial MPO (pure impurity dynamics)
	model = ImpurityHamiltonian(_rand_ham(Float64, d))
	mps4 = hybriddynamics!(sysdynamics(lattice, model, trunc=trunc), lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05))
	g0 = sysdynamics(lattice, model, trunc=trunc)
	mult!(g0, hybriddynamics(lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05)), SVDCompression(trunc))
	@test distance(mps4, g0) / norm(g0) < rtol
end


@testset "TDVPIF: PT real-time" begin

	N = 3
	δt = 0.1
	β = 1

	rtol = 5.0e-2
	trunc = truncdimcutoff(D=100, ϵ=1.0e-6, add_back=0)

	spec = Leggett(d=1, ωc=1)

	bath = bosonicbath(spec, β=β, μ=0)
	corr = Δt(bath, N=N, t=N*δt)

	d = 2
	lattice = PTLattice(N=N, δt=δt, contour=:real, d=d)

	hyb = NonAdditiveHyb(_rand_ham(ComplexF64, d))

	mps1 = hybriddynamics_naive(lattice, corr, hyb, PartialIF(trunc=trunc))

	for δ in (0.1, 0.05)
		mps2 = hybriddynamics(lattice, corr, hyb, TDVPIF(trunc=trunc, δ=δ))
		@test distance(mps1, mps2) / norm(mps1) < rtol
	end

	mps2 = hybriddynamics!(vacuumstate(ComplexF64, lattice), lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05))
	@test distance(mps1, mps2) / norm(mps1) < rtol

	# in-place flow on a non-trivial initial MPO (pure impurity dynamics)
	model = ImpurityHamiltonian(_rand_ham(ComplexF64, d))
	mps4 = hybriddynamics!(sysdynamics(lattice, model, trunc=trunc), lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05))
	g0 = sysdynamics(lattice, model, trunc=trunc)
	mult!(g0, hybriddynamics(lattice, corr, hyb, TDVPIF(trunc=trunc, δ=0.05)), SVDCompression(trunc))
	@test distance(mps4, g0) / norm(g0) < rtol
end
