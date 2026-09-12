println("------------------------------------")
println("|   Time-dependent impurity models  |")
println("------------------------------------")
println()

@testset "Time-dependent models: ADT" begin
	N = 6
	δt = 0.1
	d = 2
	T = Float64
	trunc = truncdimcutoff(D=64, ϵ=1.0e-12, add_back=0)
	tol = 1.0e-10

	ht = _rand_ham(T, d)
	hτ = _rand_ham(T, d)
	m  = _rand_ham(T, d)
	c = 2.0

	lattice_imag = ImagADTLattice(N=N, δτ=δt, d=d)
	lattice_real = RealADTLattice(N=N, δt=δt, d=d)
	lattice_mixed = MixedADTLattice(Nt=N, δt=δt, Nτ=N, δτ=δt, d=d)

	# --- quench: τ branch evolves with hτ, real branches with ht ---
	mps_ref = sysdynamics(lattice_imag, ImpurityHamiltonian(hτ))
	mps = sysdynamics(lattice_imag, QuenchedImpurityHamiltonian(hτ, ht))
	@test distance(mps, mps_ref) / norm(mps_ref) < tol

	mps_ref = sysdynamics(lattice_real, ImpurityHamiltonian(ht))
	mps = sysdynamics(lattice_real, QuenchedImpurityHamiltonian(hτ, ht))
	@test distance(mps, mps_ref) / norm(mps_ref) < tol

	for branch in (:+, :-)
		mps_ref = sysdynamics(lattice_real, ImpurityHamiltonian(ht), branch=branch)
		mps = sysdynamics(lattice_real, QuenchedImpurityHamiltonian(hτ, ht), branch=branch)
		@test distance(mps, mps_ref) / norm(mps_ref) < tol
	end

	# --- time-dependent model: τ branch evolves with the constant hτ ---
	td = TdImpurityHamiltonian(hτ, ht, [TdImpurityOp(m, t -> c * t)])
	mps_ref = sysdynamics(lattice_imag, ImpurityHamiltonian(hτ))
	mps = sysdynamics(lattice_imag, td)
	@test distance(mps, mps_ref) / norm(mps_ref) < tol

	# model(t) returns the real-time Hamiltonian ht + Σₖ fₖ(t)·mₖ
	@test isapprox(td(0.5), ht .+ (c * 0.5) .* m; rtol=0, atol=1e-14)

	# --- constant time-dependent term == the equivalent constant model ---
	td_const = TdImpurityHamiltonian(hτ, ht, [TdImpurityOp(m, t -> c)])
	mps_ref = sysdynamics(lattice_real, ImpurityHamiltonian(ht .+ c .* m))
	mps = sysdynamics(lattice_real, td_const)
	@test distance(mps, mps_ref) / norm(mps_ref) < tol

	# --- genuinely time-dependent: per-step propagators evaluated at
	#     t = (j-1)δt (forward) resp. t = (N-j)δt (backward); compared
	#     against a manual sequential application of the exact propagators ---
	function manual_td_reference(lattice, branch, hτ_, ht_, m_, f, N, trunc)
		g = vacuumstate(ComplexF64, lattice)
		alg = Orthogonalize(SVD(), trunc)
		for j in 1:N
			t = (branch == :+) ? (j - 1) * lattice.δt : (N - j) * lattice.δt
			U = exp((branch == :+ ? -im : im) * lattice.δt .* (ht_ .+ f(t) .* m_))
			a, b = (branch == :-) ? (j, j+1) : (j+1, j)
			term = ADTTerm((index(lattice, a, branch=branch), index(lattice, b, branch=branch)), U)
			apply!(term, g)
			canonicalize!(g, alg=alg)
		end
		return g
	end

	f = t -> 0.5 + 0.25 * t
	td = TdImpurityHamiltonian(hτ, ht, [TdImpurityOp(m, f)])
	for branch in (:+, :-)
		mps_ref = manual_td_reference(lattice_real, branch, hτ, ht, m, f, N, trunc)
		mps = sysdynamics(lattice_real, td, branch=branch)
		@test distance(mps, mps_ref) / norm(mps_ref) < tol
	end

	# --- mixed contour routing: constant Td term == equivalent quench model on all branches ---
	mps_ref = sysdynamics(lattice_mixed, QuenchedImpurityHamiltonian(hτ, ht .+ c .* m))
	mps = sysdynamics(lattice_mixed, td_const)
	@test distance(mps, mps_ref) / norm(mps_ref) < tol
end

@testset "Time-dependent models: PT" begin
	N = 6
	δt = 0.1
	d = 2
	T = Float64
	trunc = truncdimcutoff(D=64, ϵ=1.0e-12, add_back=0)
	tol = 1.0e-10

	ht = _rand_ham(T, d)
	hτ = _rand_ham(T, d)
	m  = _rand_ham(T, d)
	c = 2.0

	lattice_imag = ImagPTLattice(N=N, δτ=δt, d=d)
	lattice_real = RealPTLattice(N=N, δt=δt, d=d)
	lattice_mixed = MixedPTLattice(Nt=N, δt=δt, Nτ=N, δτ=δt, d=d)

	# --- quench ---
	mps_ref = sysdynamics(lattice_imag, ImpurityHamiltonian(hτ))
	mps = sysdynamics(lattice_imag, QuenchedImpurityHamiltonian(hτ, ht))
	@test distance(mps, mps_ref) / norm(mps_ref) < tol

	mps_ref = sysdynamics(lattice_real, ImpurityHamiltonian(ht))
	mps = sysdynamics(lattice_real, QuenchedImpurityHamiltonian(hτ, ht))
	@test distance(mps, mps_ref) / norm(mps_ref) < tol

	for branch in (:+, :-)
		mps_ref = sysdynamics(lattice_real, ImpurityHamiltonian(ht), branch=branch)
		mps = sysdynamics(lattice_real, QuenchedImpurityHamiltonian(hτ, ht), branch=branch)
		@test distance(mps, mps_ref) / norm(mps_ref) < tol
	end

	# --- time-dependent model ---
	td_const = TdImpurityHamiltonian(hτ, ht, [TdImpurityOp(m, t -> c)])
	mps_ref = sysdynamics(lattice_imag, ImpurityHamiltonian(hτ))
	mps = sysdynamics(lattice_imag, td_const)
	@test distance(mps, mps_ref) / norm(mps_ref) < tol

	# constant time-dependent term == the equivalent constant model
	mps_ref = sysdynamics(lattice_real, ImpurityHamiltonian(ht .+ c .* m))
	mps = sysdynamics(lattice_real, td_const)
	@test distance(mps, mps_ref) / norm(mps_ref) < tol

	# genuinely time-dependent, manual per-step reference on the forward branch
	f = t -> 0.5 + 0.25 * t
	td = TdImpurityHamiltonian(hτ, ht, [TdImpurityOp(m, f)])
	g = vacuumstate(ComplexF64, lattice_real)
	alg = Orthogonalize(SVD(), trunc)
	for j in 1:N
		t = (j - 1) * lattice_real.δt
		U = exp(-im * lattice_real.δt .* (ht .+ f(t) .* m))
		apply!(ContourOperator(ContourIndex(j, :+), U), lattice_real, g)
	end
	canonicalize!(g, alg=alg)
	mps = sysdynamics(lattice_real, td, branch=:+)
	@test distance(mps, g) / norm(g) < tol

	# mixed contour routing
	mps_ref = sysdynamics(lattice_mixed, QuenchedImpurityHamiltonian(hτ, ht .+ c .* m))
	mps = sysdynamics(lattice_mixed, td_const)
	@test distance(mps, mps_ref) / norm(mps_ref) < tol
end