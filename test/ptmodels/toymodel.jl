println("------------------------------------")
println("|              Toy Model           |")
println("------------------------------------")

# H = Ω*σz	+ σx(a + a†) + σy(b + b†) + a† a + b† b
@testset "Toy model: imaginary-time" begin

	Ω = 0.5
	N = 10
	δτ = 0.05
	β = N * δτ
	chi = 100
	d = 20
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(N = N, δτ=δτ, contour=:imag)

	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]
	hop = Ω .* z
	Is = one(x)
	Ib = one(zeros(d, d))
	model = ImpurityHamiltonian(hop)

	nop = bosondensityoperator(d=d)
	Hbarebath = kron(nop, Ib) + kron(Ib, nop)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib, Ib) + kron(Is, Hbarebath) + kron(x, a' + a, Ib) + kron(y, Ib, a' + a)

	mpsK = sysdynamics(lattice, model, trunc=trunc)
	
	bs1 = NonAdditiveHyb(x)
	bs2 = NonAdditiveHyb(y)

	spec = DiracDelta(1)

	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)

	# mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI1 = hybriddynamics_naive(lattice, corr, bs1, trunc=trunc)
	mpsI2 = hybriddynamics_naive(lattice, corr, bs2, trunc=trunc)
	# println(distance(mpsI, mpsI′), " ", norm(mpsI), " ", norm(mpsI′))
	# @test distance(mpsI, mpsI′) / norm(mpsI′) < tol
	mps = mult!(complex(mpsK), mpsI1, trunc=trunc)
	mps = mult!(mps, mpsI2, trunc=trunc)

	Zval = integrate(lattice, mps)

	ρ = exp(-β * H)
	Zval2 = tr(ρ) / tr(exp(-β .* Hbarebath))
	@test abs(Zval - Zval2) / abs(Zval) < tol


	## diagonal observables
	op = [-0.73 0; 0 0.5]

	ind1 = ContourIndex(1)
	t = ContourOperator(ind1, op * op )
	mps2 = apply!(t, lattice, deepcopy(mps))
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	for i in 2:N
		ind2 = ContourIndex(i)
		t = ContourOperator([ind2,ind1], [op, op])
		# t = ADTTerm((i,1), reshape(kron(zdiag, zdiag), 2, 2))
		mps2 = apply!(t, lattice, deepcopy(mps))
		v = integrate(lattice, mps2) / Zval
		push!(corrs, v)
	end
	

	A = kron(op, Ib, Ib)

	corrs2 = correlation_2op_1τ(H, A, A, 0:δτ:β, β=β)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	## off-diagonal observables
	op1 = [0 0; 0.7 0]
	op2 = [0 0.8;0 0 ]

	c1 = ContourIndex(1)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = apply!(ct, lattice, deepcopy(mps))
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	c2 = ContourIndex(1)
	for i in 2:N
		c1 = ContourIndex(i)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = apply!(ct, lattice, deepcopy(mps))
		v = integrate(lattice, mps2) / Zval
		push!(corrs, v)
	end
	

	A1 = kron(op1, Ib, Ib)
	A2 = kron(op2, Ib, Ib)

	corrs2 = correlation_2op_1τ(H, A1, A2, 0:δτ:β, β=β)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end


@testset "Toy model: real-time" begin

	Ω = 0.5
	N = 10
	δt = 0.05
	β = 2
	t = N * δt
	chi = 100
	d = 20
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(N = N, δt=δt, contour=:real)

	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]
	hop = Ω .* z
	Is = one(x)
	Ib = one(zeros(d, d))
	model = ImpurityHamiltonian(hop)

	nop = bosondensityoperator(d=d)
	Hbarebath = kron(nop, Ib) + kron(Ib, nop)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib, Ib) + kron(Is, Hbarebath) + kron(x, a' + a, Ib) + kron(y, Ib, a' + a)


	bs1 = NonAdditiveHyb(x)
	bs2 = NonAdditiveHyb(y)


	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	# mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI1 = hybriddynamics_naive(lattice, corr, bs1, trunc=trunc)
	mpsI2 = hybriddynamics_naive(lattice, corr, bs2, trunc=trunc)

	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mps = mult!(mpsK, mpsI1, trunc=trunc)
	mps = mult!(mps, mpsI2, trunc=trunc)

	ρ1 = zeros(2,2)
	ρ1[1,1] = 1
	ρ2 = 0.5 .* one(ρ1)

	for ρimp in [ρ1, _rand_dm(2)]

		tmp = initialstate!(deepcopy(mps), lattice, ρimp)
		Zval = integrate(lattice, tmp)
		

		
		ρ = kron(ρimp, exp(-β * Hbarebath)) 

		## diagonal observables
		op = [-0.73 0; 0 0.5]

		ind1 = ContourIndex(1, branch=:+)
		m = ContourOperator(ind1, op * op )
		mps2 = apply!(m, lattice, deepcopy(mps))
		mps2 = initialstate!(mps2, lattice, ρimp)
		v = integrate(lattice, mps2) / Zval

		corrs = [v]
		for i in 2:N
			ind2 = ContourIndex(i, branch=:+)
			m = ContourOperator([ind2,ind1], [op, op])
			mps2 = apply!(m, lattice, deepcopy(mps))
			mps2 = initialstate!(mps2, lattice, ρimp)
			v = integrate(lattice, mps2) / Zval
			push!(corrs, v)
		end
		
		A = kron(op, Ib, Ib)
		corrs2 = correlation_2op_1t(H, A, A, ρ, 0:δt:t, reverse = false)
		corrs2 = corrs2[1:length(corrs)]
		@test norm(corrs - corrs2) / norm(corrs2) < tol


		# off-diagonal observables

		op1 = [0 0.8; 0 0]
		op2 = [0 0; 0.7*im 0]

		A1 = kron(op1, Ib, Ib)
		A2 = kron(op2, Ib, Ib)


		c1 = ContourIndex(1, branch=:+)

		ct = ContourOperator(c1, op1 * op2)
		mps2 = apply!(ct, lattice, deepcopy(mps))
		mps2 = initialstate!(mps2, lattice, ρimp)
		v = integrate(lattice, mps2) / Zval

		corrs = [v]
		c2 = ContourIndex(1, branch=:+)
		for i in 2:N
			c1 = ContourIndex(i, branch=:+)
			ct = ContourOperator([c1, c2], [op1, op2])

			mps2 = apply!(ct, lattice, deepcopy(mps))
			mps2 = initialstate!(mps2, lattice, ρimp)
			v = integrate(lattice, mps2) / Zval

			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = false)
		corrs2 = corrs2[1:length(corrs)]

		@test norm(corrs - corrs2) / norm(corrs2) < tol

		op1 = [0 0.8*im; 0 0]
		op2 = [0 0; 0.7 0]


		A1 = kron(op1, Ib, Ib)
		A2 = kron(op2, Ib, Ib)

		c1 = ContourIndex(1, branch=:-)

		ct = ContourOperator(c1, op1 * op2)
		mps2 = apply!(ct, lattice, deepcopy(mps))
		mps2 = initialstate!(mps2, lattice, ρimp)
		v = integrate(lattice, mps2) / Zval

		corrs = [v]
		for i in 2:N
			c2 = ContourIndex(i, branch=:+)
			ct = ContourOperator([c1, c2], [op1, op2])

			mps2 = apply!(ct, lattice, deepcopy(mps))
			mps2 = initialstate!(mps2, lattice, ρimp)
			v = integrate(lattice, mps2) / Zval

			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = true)
		corrs2 = corrs2[1:length(corrs)]

		@test norm(corrs - corrs2) / norm(corrs2) < tol

	end

end

@testset "Rabi model: mixed-time" begin

	Ω = 0.5
	Nt = 5
	δt = 0.03
	t = Nt * δt
	Nτ = 10
	δτ = 0.05
	β = Nτ * δτ
	chi = 100
	d = 20
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(Nt = Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed)

	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]
	hop = Ω .* z
	Is = one(x)
	Ib = one(zeros(d, d))
	model = ImpurityHamiltonian(hop)

	nop = bosondensityoperator(d=d)
	Hbarebath = kron(nop, Ib) + kron(Ib, nop)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib, Ib) + kron(Is, Hbarebath) + kron(x, a' + a, Ib) + kron(y, Ib, a' + a)


	op1 = [0 0.8; 0 0]
	op2 = [0 0; 0.7 0]


	A1 = kron(op1, Ib, Ib)
	A2 = kron(op2, Ib, Ib)


	bs1 = NonAdditiveHyb(x)
	bs2 = NonAdditiveHyb(y)


	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	# mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI1 = hybriddynamics_naive(lattice, corr, bs1, trunc=trunc)
	mpsI2 = hybriddynamics_naive(lattice, corr, bs2, trunc=trunc)
	# @test distance(mpsI, mpsI′) / norm(mpsI′) < tol

	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mps = mult!(mpsK, mpsI1, trunc=trunc)
	mps = mult!(mps, mpsI2, trunc=trunc)
	Zval = integrate(lattice, mps)


	# H, Hbarebath = rabi_ham_2(Ω, d=d)


	ρ = exp(-β .* H)

	# off-diagonal observables

	c1 = ContourIndex(1, branch=:+)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = apply!(ct, lattice, deepcopy(mps))
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	c2 = ContourIndex(1, branch=:+)
	for i in 2:Nt
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = apply!(ct, lattice, deepcopy(mps))
		v = integrate(lattice, mps2) / Zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = false)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	c1 = ContourIndex(1, branch=:-)

	ct = ContourOperator(c1, op2 * op1)
	mps2 = apply!(ct, lattice, deepcopy(mps))
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	for i in 2:Nt
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op2, op1])

		mps2 = apply!(ct, lattice, deepcopy(mps))
		v = integrate(lattice, mps2) / Zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(H, A2, A1, ρ, 0:δt:t, reverse = true)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end

