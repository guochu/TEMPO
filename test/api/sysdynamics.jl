# Impurity dynamics on ADT / PT lattices: boundary conditions, sysdynamics! variants,
# quenched and explicitly time-dependent models

@testset "boundarycondition          " begin
	lattice = ADTLattice(N=3, δt=0.1, contour=:real)
	x = randomadt(ComplexF64, length(lattice), D=4, d=2)
	canonicalize!(x)
	xc = copy(x)
	ρimp = spin_up()
	x2 = boundarycondition(x, lattice, ρ₀=ρimp)
	@test x2 !== x
	@test distance(x2, boundarycondition!(copy(x), lattice, ρ₀=ρimp)) < 1.0e-12
	# the non-mutating version must not modify the input
	@test distance(x, xc) < 1.0e-14
end

@testset "sysdynamics! branches      " begin
	N = 4
	δt = 0.1
	lattice = ADTLattice(N=N, δt=δt, contour=:real)
	model = ImpurityHamiltonian(0.5 .* pauli_x())
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)
	a = sysdynamics!(vacuumstate(ComplexF64, lattice), lattice, model, trunc=trunc)
	# forward branch only differs from both-branch evolution
	b = sysdynamics!(vacuumstate(ComplexF64, lattice), lattice, model, trunc=trunc, branch=:+)
	@test distance(a, b) > 1.0e-6
	sysdynamics!(b, lattice, model, trunc=trunc, branch=:-)
	@test distance(a, b) < 1.0e-10
	# non-mutating wrapper consistency
	m1 = sysdynamics(lattice, model, trunc=trunc)
	m2 = sysdynamics!(vacuumstate(ComplexF64, lattice), lattice, model, trunc=trunc)
	@test distance(m1, m2) < 1.0e-12
end

@testset "Time-dependent impurity models" begin
	N = 6
	δt = 0.1
	d = 2
	trunc = truncdimcutoff(D=64, ϵ=1.0e-12, add_back=0)
	tol = 1.0e-10

	for name in ("ADT", "PT")
		ht = _rand_ham(Float64, d)
		hτ = _rand_ham(Float64, d)
		m  = _rand_ham(Float64, d)
		c = 2.0

		lattice_imag, lattice_real, lattice_mixed = if name == "ADT"
			(ImagADTLattice(N=N, δτ=δt, d=d), RealADTLattice(N=N, δt=δt, d=d), MixedADTLattice(Nt=N, δt=δt, Nτ=N, δτ=δt, d=d))
		else
			(ImagPTLattice(N=N, δτ=δt, d=d), RealPTLattice(N=N, δt=δt, d=d), MixedPTLattice(Nt=N, δt=δt, Nτ=N, δτ=δt, d=d))
		end

		@testset "$name" begin
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

			# --- genuinely time-dependent vs a manual sequential application of
			#     the exact per-step propagators on the forward branch ---
			f = t -> 0.5 + 0.25 * t
			td = TdImpurityHamiltonian(hτ, ht, [TdImpurityOp(m, f)])
			alg = Orthogonalize(SVD(), trunc)
			g = vacuumstate(ComplexF64, lattice_real)
			for j in 1:N
				t = (j - 1) * lattice_real.δt
				U = exp(-im * lattice_real.δt .* (ht .+ f(t) .* m))
				if name == "ADT"
					term = ADTTerm((index(lattice_real, j + 1, branch=:+), index(lattice_real, j, branch=:+)), U)
					apply!(term, g)
					canonicalize!(g, alg=alg)
				else
					apply!(ContourOperator(ContourIndex(j, :+), U), lattice_real, g)
				end
			end
			name == "ADT" || canonicalize!(g, alg=alg)
			mps = sysdynamics(lattice_real, td, branch=:+)
			@test distance(mps, g) / norm(g) < tol

			# --- mixed contour routing: constant Td term == equivalent quench model on all branches ---
			mps_ref = sysdynamics(lattice_mixed, QuenchedImpurityHamiltonian(hτ, ht .+ c .* m))
			mps = sysdynamics(lattice_mixed, td_const)
			@test distance(mps, mps_ref) / norm(mps_ref) < tol
		end
	end
end
