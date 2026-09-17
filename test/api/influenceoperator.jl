# influenceoperators / influenceoperatorsteppers / influenceoperatorstepper vs manual
# accumulation of the influence operator, parametrized over ADT (ADTTerm) and PT (FockTermS)

@testset "InfluenceOperator: ADT lattice" begin
	# --- imaginary time ---
	β = 1
	δτ = 0.2
	tol = 1.0e-6
	N = round(Int, β/δτ)
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=50, ϵ=1.0e-6, add_back=0)
	algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
	d = 2

	lattice = ADTLattice(N=N, δτ=δτ, d=d, contour=:imag)
	corr = correlationfunction(bath, lattice)
	hyb = AdditiveHyb(randn(Float64, d))

	z = hyb.op
	z2 = z .* z
	zz = reshape(kron(z, z), d, d)

	mpo1 = only(influenceoperators(lattice, corr, hyb, algexpan=algexpan))

	orth = Orthogonalize(SVD(), trunc)
	mpo2 = nothing
	for i in 1:lattice.N, j in 1:lattice.N
		ind1, ind2 = ContourIndex(i+1), ContourIndex(j+1)
		coef = index(corr, i, j)
		if coef != 0
			if ind1 == ind2
				t = ADTTerm(lattice[ind1], coef .* z2)
			else
				t = ADTTerm((lattice[ind1], lattice[ind2]), coef .* zz)
			end
			if isnothing(mpo2)
				mpo2 = apply!(t, vacuumstate(Float64, lattice))
			else
				mpo2 += apply!(t, vacuumstate(Float64, lattice))
			end
			canonicalize!(mpo2, alg=orth)
		end
	end
	@test distance(mpo1, mpo2) / norm(mpo2) < tol

	# one first-order step: I + dt·K
	dt = 0.01
	mpo1 = dt * mpo1 + vacuumstate(Float64, lattice)
	mps0, = influenceoperatorsteppers(lattice, corr, dt, hyb, WII(), algexpan=algexpan)
	@test distance(mpo1, mps0) / norm(mps0) < dt

	# single steppers for all four stepper types and two multiplication algorithms
	for algmult in (SVDCompression(truncdimcutoff(D=50, ϵ=1.0e-12, add_back=0)), DMRG1(trunc=truncdimcutoff(D=50, ϵ=1.0e-6)))
		mps1 = influenceoperatorstepper(lattice, corr, dt, hyb, WII(), algmult, algexpan=algexpan)
		_n = norm(mps1)
		mps2 = influenceoperatorstepper(lattice, corr, dt, hyb, WI(), algmult, algexpan=algexpan)
		mps3 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WI()), algmult, algexpan=algexpan)
		mps4 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WII()), algmult, algexpan=algexpan)
		@test distance(mps1, mps0) / _n < dt
		@test distance(mps2, mps0) / _n < dt
		@test distance(mps3, mps0) / _n < dt
		@test distance(mps4, mps0) / _n < dt
	end

	# --- real time: one partial influence operator per branch combination ---
	β = 1.0
	δt = 0.1
	N = 4
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-8, add_back=0)
	algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
	tol = 1.0e-5

	lattice = ADTLattice(N=N, δt=δt, contour=:real)
	corr = correlationfunction(bath, lattice)
	hyb = AdditiveHyb(randn(Float64, d))
	z = hyb.op
	z2 = z .* z
	zz = reshape(kron(z, z), d, d)

	function manual_if(b1, b2)
		mps = nothing
		for i in 1:lattice.N, j in 1:lattice.N
			ind1, ind2 = ContourIndex(i, b1), ContourIndex(j, b2)
			coef = index(corr, i, j, b1=b1, b2=b2)
			if coef != 0
				if ind1 == ind2
					t = ADTTerm(lattice[ind1], coef .* z2)
				else
					t = ADTTerm((lattice[ind1], lattice[ind2]), coef .* zz)
				end
				if isnothing(mps)
					mps = apply!(t, vacuumstate(ComplexF64, lattice))
				else
					mps += apply!(t, vacuumstate(ComplexF64, lattice))
				end
			end
		end
		return mps
	end

	ifs = influenceoperators(lattice, corr, hyb, algexpan=algexpan)
	@test length(ifs) == 4
	for (b1, b2) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
		mps_ref = manual_if(b1, b2)
		mps = (b1 == :+ && b2 == :+) ? ifs[1] : (b1 == :+ ? ifs[2] : (b2 == :+ ? ifs[3] : ifs[4]))
		@test distance(mps, mps_ref) / norm(mps_ref) < tol
	end

	# steppers: one first-order step for each branch combination
	dt = 0.01
	mps_ref = [dt * m + vacuumstate(ComplexF64, lattice) for m in ifs]
	mps_stp = influenceoperatorsteppers(lattice, corr, dt, hyb, WII(), algexpan=algexpan)
	for (m1, m2) in zip(mps_ref, mps_stp)
		@test distance(m1, m2) / norm(m2) < dt
	end

	mps0 = mult!(mps_stp[1], mps_stp[2], trunc=trunc)
	mps0 = mult!(mps0, mps_stp[3], trunc=trunc)
	mps0 = mult!(mps0, mps_stp[4], trunc=trunc)
	for algmult in (SVDCompression(truncdimcutoff(D=50, ϵ=1.0e-12, add_back=0)), DMRG1(trunc=truncdimcutoff(D=50, ϵ=1.0e-6)))
		mps1 = influenceoperatorstepper(lattice, corr, dt, hyb, WII(), algmult, algexpan=algexpan)
		_n = norm(mps1)
		mps2 = influenceoperatorstepper(lattice, corr, dt, hyb, WI(), algmult, algexpan=algexpan)
		mps3 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WI()), algmult, algexpan=algexpan)
		mps4 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WII()), algmult, algexpan=algexpan)
		@test distance(mps1, mps0) / _n < dt
		@test distance(mps2, mps0) / _n < dt
		@test distance(mps3, mps0) / _n < dt
		@test distance(mps4, mps0) / _n < dt
	end
end

@testset "InfluenceOperator: PT lattice" begin
	__ti_localop(η, op1, op2) = η.ηⱼₖ[1] .* (op1 * op2) + η.ηₖⱼ[1] .* (op2 * op1)

	# --- imaginary time ---
	β = 1
	δτ = 0.2
	tol = 1.0e-6
	N = round(Int, β/δτ)
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=50, ϵ=1.0e-6, add_back=0)
	algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
	d = 2

	lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)
	corr = correlationfunction(bath, lattice)

	for hyb in (NonAdditiveHyb(_rand_ham(d)), NonDiagonalHyb(randn(ComplexF64, d, d)))
		op1, op2 = pairop(hyb)

		mpo1 = only(influenceoperators(lattice, corr, hyb, algexpan=algexpan))

		mpo2 = nothing
		for i in 1:lattice.N, j in 1:lattice.N
			ind1, ind2 = ContourIndex(i), ContourIndex(j)
			if ind1 == ind2
				m = __ti_localop(corr.data, op1, op2)
				t = FockTermS(lattice[ind1], m)
			else
				coef = index(corr, i, j)
				m = coef .* kron(op2, op1)
				t = FockTermS((lattice[ind1], lattice[ind2]), reshape(m, (d, d, d, d)))
			end
			if isnothing(mpo2)
				mpo2 = apply!(t, vacuumstate(scalartype(hyb), lattice))
			else
				mpo2 += apply!(t, vacuumstate(scalartype(hyb), lattice))
			end
		end
		canonicalize!(mpo2, alg=Orthogonalize(trunc=trunc))
		@test distance(mpo1, mpo2) / norm(mpo2) < tol

		dt = 0.01
		mpo1 = dt * mpo1 + vacuumstate(Float64, lattice)
		mps0, = influenceoperatorsteppers(lattice, corr, dt, hyb, WII(), algexpan=algexpan)
		@test distance(mpo1, mps0) / norm(mps0) < dt

		for algmult in (SVDCompression(truncdimcutoff(D=50, ϵ=1.0e-12, add_back=0)), DMRG1(trunc=truncdimcutoff(D=50, ϵ=1.0e-6)))
			mps1 = influenceoperatorstepper(lattice, corr, dt, hyb, WII(), algmult, algexpan=algexpan)
			_n = norm(mps1)
			mps2 = influenceoperatorstepper(lattice, corr, dt, hyb, WI(), algmult, algexpan=algexpan)
			mps3 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WI()), algmult, algexpan=algexpan)
			mps4 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WII()), algmult, algexpan=algexpan)
			@test distance(mps1, mps0) / _n < dt
			@test distance(mps2, mps0) / _n < dt
			@test distance(mps3, mps0) / _n < dt
			@test distance(mps4, mps0) / _n < dt
		end
	end

	# --- real time ---
	__get_contour_op(lattice, ind1, ind2, z1, z2, corr) = begin
		d = lattice.d
		(branch(ind1) == :-) && (z1 = transpose(z1))
		(branch(ind2) == :-) && (z2 = transpose(z2))
		η = branch(corr, branch(ind1), branch(ind2))
		pos1, pos2 = lattice[ind1], lattice[ind2]
		if ind1 == ind2
			m = __ti_localop(η, z1, z2)
			t = FockTermS(pos1, m)
		else
			coef = η[ind1.j, ind2.j]
			m = coef .* kron(z2, z1)
			t = FockTermS((pos1, pos2), reshape(m, (d, d, d, d)))
		end
		return t
	end

	β = 1.0
	δt = 0.1
	N = 4
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-8, add_back=0)
	algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
	tol = 1.0e-5

	lattice = PTLattice(N=N, δt=δt, contour=:real)
	corr = correlationfunction(bath, lattice)

	for hyb in (NonAdditiveHyb(_rand_ham(d)), NonDiagonalHyb(randn(ComplexF64, d, d)))
		op1, op2 = pairop(hyb)

		function manual_if(b1, b2)
			mps = nothing
			for i in 1:lattice.N, j in 1:lattice.N
				ind1, ind2 = ContourIndex(i, b1), ContourIndex(j, b2)
				t = __get_contour_op(lattice, ind1, ind2, op1, op2, corr)
				if isnothing(mps)
					mps = apply!(t, vacuumstate(scalartype(hyb), lattice))
				else
					mps += apply!(t, vacuumstate(scalartype(hyb), lattice))
				end
			end
			return mps
		end

		ifs = influenceoperators(lattice, corr, hyb, algexpan=algexpan)
		@test length(ifs) == 4
		for (k, (b1, b2)) in enumerate(((:+, :+), (:+, :-), (:-, :+), (:-, :-)))
			mps_ref = manual_if(b1, b2)
			@test distance(ifs[k], mps_ref) / norm(mps_ref) < tol
		end

		dt = 0.01
		mps_ref = [dt * m + vacuumstate(ComplexF64, lattice) for m in ifs]
		mps_stp = influenceoperatorsteppers(lattice, corr, dt, hyb, WII(), algexpan=algexpan)
		for (m1, m2) in zip(mps_ref, mps_stp)
			@test distance(m1, m2) / norm(m2) < dt
		end

		mps0 = mult!(mps_stp[1], mps_stp[2], trunc=trunc)
		mps0 = mult!(mps0, mps_stp[3], trunc=trunc)
		mps0 = mult!(mps0, mps_stp[4], trunc=trunc)
		for algmult in (SVDCompression(truncdimcutoff(D=50, ϵ=1.0e-12, add_back=0)), DMRG1(trunc=truncdimcutoff(D=50, ϵ=1.0e-6)))
			mps1 = influenceoperatorstepper(lattice, corr, dt, hyb, WII(), algmult, algexpan=algexpan)
			_n = norm(mps1)
			mps2 = influenceoperatorstepper(lattice, corr, dt, hyb, WI(), algmult, algexpan=algexpan)
			mps3 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WI()), algmult, algexpan=algexpan)
			mps4 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WII()), algmult, algexpan=algexpan)
			@test distance(mps1, mps0) / _n < dt
			@test distance(mps2, mps0) / _n < dt
			@test distance(mps3, mps0) / _n < dt
			@test distance(mps4, mps0) / _n < dt
		end
	end
end
