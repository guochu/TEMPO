println("------------------------------------")
println("|            Rabi Model            |")
println("------------------------------------")


@testset "Rabi model: imaginary-time" begin

	Ω = 0.5
	N = 20
	δτ = 0.1
	β = N * δτ
	chi = 100
	d = 50
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(N = N, δτ=δτ, contour=:imag)

	x = [0 1; 1 0]
	hop = Ω .* x
	z = [-1 0; 0 1]
	Is = one(x)
	Ib = one(zeros(d, d))
	model = ImpurityHamiltonian(hop)

	mpsK = sysdynamics(lattice, model, trunc=trunc)
	
	bs = NonAdditiveHyb(z)

	spec = DiracDelta(1)

	bath = bosonicbath(spec, β=β)

	corr = correlationfunction(bath, lattice)

	# mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI = hybriddynamics_naive(lattice, corr, bs, trunc=trunc)
	# println(distance(mpsI, mpsI′), " ", norm(mpsI), " ", norm(mpsI′))
	# @test distance(mpsI, mpsI′) / norm(mpsI′) < tol
	mps = mult!(mpsK, mpsI, trunc=trunc)


	H, Hbarebath = rabi_ham(Ω, d=d)

	ρ = exp(-β * H)
	Zval = integrate(lattice, mps)
	@test abs(Zval - tr(ρ) / tr(exp(-β .* Hbarebath))) / abs(Zval) < tol


	## diagonal observables
	op = [-0.73 0; 0 0.5]

	ind1 = ContourIndex(1)
	t = ContourOperator(ind1, op * op )
	mps2 = apply!(t, lattice, deepcopy(mps))
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	ids2N = [k for k in sampleidx(N) if k > 1]
	idsallN = [1; ids2N]
	for i in ids2N
		ind2 = ContourIndex(i)
		t = ContourOperator([ind2,ind1], [op, op])
		# t = ADTTerm((i,1), reshape(kron(zdiag, zdiag), 2, 2))
		mps2 = apply!(t, lattice, deepcopy(mps))
		v = integrate(lattice, mps2) / Zval
		push!(corrs, v)
	end
	

	A = kron(op, Ib)

	corrs2 = correlation_2op_1τ(H, A, A, 0:δτ:β, β=β)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	## off-diagonal observables
	op1 = [0 0; 0.7 0]
	op2 = [0 0.8;0 0 ]

	c1 = ContourIndex(1)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = apply!(ct, lattice, deepcopy(mps))
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	for i in ids2N
		c2 = ContourIndex(i)
		ct = ContourOperator([c2, c1], [op2, op1])

		mps2 = apply!(ct, lattice, deepcopy(mps))
		v = integrate(lattice, mps2) / Zval
		push!(corrs, v)
	end
	

	A1 = kron(op1, Ib)
	A2 = kron(op2, Ib)

	corrs2 = correlation_2op_1τ(H, A1, A2, 0:δτ:β, β=β)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end

@testset "Rabi model: real-time" begin

	Ω = 0.5
	N = 10
	δt = 0.05
	β = 2
	t = N * δt
	chi = 100
	d = 50
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(N = N, δt=δt, contour=:real)

	# x = [0 1; 1 0]
	x = Matrix{ComplexF64}([0 im; -im 0])
	hop = Ω .* x
	z = [-1 0; 0 1]
	Is = one(x)
	Ib = one(zeros(d, d))
	model = ImpurityHamiltonian(hop)

	Hbarebath = bosondensityoperator(d=d)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib) + kron(Is, Hbarebath) + kron(z, a' + a)


	bs = NonAdditiveHyb(z)
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	# mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI = hybriddynamics_naive(lattice, corr, bs, trunc=trunc)
	# @test distance(mpsI, mpsI′) / norm(mpsI′) < tol
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mps = mult!(mpsK, mpsI, trunc=trunc)

	ρ1 = zeros(2,2)
	ρ1[1,1] = 1
	ρ2 = 0.5 .* one(ρ1)

	for ρimp in [ρ1, ρ2]

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
		ids2N = [k for k in sampleidx(N) if k > 1]
		idsallN = [1; ids2N]
		for i in ids2N
			ind2 = ContourIndex(i, branch=:+)
			m = ContourOperator([ind2,ind1], [op, op])
			mps2 = apply!(m, lattice, deepcopy(mps))
			mps2 = initialstate!(mps2, lattice, ρimp)
			v = integrate(lattice, mps2) / Zval
			push!(corrs, v)
		end
		
		A = kron(op, Ib)
		corrs2 = correlation_2op_1t(H, A, A, ρ, 0:δt:t, reverse = false)
		corrs2 = corrs2[idsallN]
		@test norm(corrs - corrs2) / norm(corrs2) < tol


		# off-diagonal observables

		op1 = [0 0.8; 0 0]
		op2 = [0 0; 0.7*im 0]

		A1 = kron(op1, Ib)
		A2 = kron(op2, Ib)


		c1 = ContourIndex(1, branch=:+)

		ct = ContourOperator(c1, op1 * op2)
		mps2 = apply!(ct, lattice, deepcopy(mps))
		mps2 = initialstate!(mps2, lattice, ρimp)
		v = integrate(lattice, mps2) / Zval

		corrs = [v]
		c2 = ContourIndex(1, branch=:+)
		for i in ids2N
			c1 = ContourIndex(i, branch=:+)
			ct = ContourOperator([c1, c2], [op1, op2])

			mps2 = apply!(ct, lattice, deepcopy(mps))
			mps2 = initialstate!(mps2, lattice, ρimp)
			v = integrate(lattice, mps2) / Zval

			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = false)
		corrs2 = corrs2[idsallN]

		@test norm(corrs - corrs2) / norm(corrs2) < tol

		op1 = [0 0.8*im; 0 0]
		op2 = [0 0; 0.7 0]


		A1 = kron(op1, Ib)
		A2 = kron(op2, Ib)

		c1 = ContourIndex(1, branch=:-)

		ct = ContourOperator(c1, op1 * op2)
		mps2 = apply!(ct, lattice, deepcopy(mps))
		mps2 = initialstate!(mps2, lattice, ρimp)
		v = integrate(lattice, mps2) / Zval

		corrs = [v]
		for i in ids2N
			c2 = ContourIndex(i, branch=:+)
			ct = ContourOperator([c1, c2], [op1, op2])

			mps2 = apply!(ct, lattice, deepcopy(mps))
			mps2 = initialstate!(mps2, lattice, ρimp)
			v = integrate(lattice, mps2) / Zval

			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = true)
		corrs2 = corrs2[idsallN]

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
	d = 50
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(Nt = Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed)

	x = [0 1; 1 0]
	hop = Ω .* x
	z = [-1 0; 0 1]
	Is = one(x)
	Ib = one(zeros(d, d))
	model = ImpurityHamiltonian(hop)

	op1 = [0 0.8; 0 0]
	op2 = [0 0; 0.7 0]


	A1 = kron(op1, Ib)
	A2 = kron(op2, Ib)


	bs = NonAdditiveHyb(z)
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	# mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI = hybriddynamics_naive(lattice, corr, bs, trunc=trunc)
	# @test distance(mpsI, mpsI′) / norm(mpsI′) < tol

	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mps = mult!(mpsK, mpsI, trunc=trunc)
	Zval = integrate(lattice, mps)


	H, Hbarebath = rabi_ham(Ω, d=d)


	ρ = exp(-β .* H)

	# off-diagonal observables

	c1 = ContourIndex(1, branch=:+)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = apply!(ct, lattice, deepcopy(mps))
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	c2 = ContourIndex(1, branch=:+)
	ids2Nt = [k for k in sampleidx(Nt) if k > 1]
	idsallNt = [1; ids2Nt]
	for i in ids2Nt
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = apply!(ct, lattice, deepcopy(mps))
		v = integrate(lattice, mps2) / Zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = false)
	corrs2 = corrs2[idsallNt]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	c1 = ContourIndex(1, branch=:-)

	ct = ContourOperator(c1, op2 * op1)
	mps2 = apply!(ct, lattice, deepcopy(mps))
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	for i in ids2Nt
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op2, op1])

		mps2 = apply!(ct, lattice, deepcopy(mps))
		v = integrate(lattice, mps2) / Zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(H, A2, A1, ρ, 0:δt:t, reverse = true)
	corrs2 = corrs2[idsallNt]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end
println("------------------------------------")
println("|     Dissipative Rabi Model       |")
println("------------------------------------")

@testset "Rabi model: real-time (dissipative)" begin

	Ω = 0.5
	N = 6  # Liouville-space dynamics (d² = 4): keep the chain short
	δt = 0.05
	β = 2
	t = N * δt
	chi = 100
	d = 20
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(N = N, δt=δt, contour=:real)

	# x = [0 1; 1 0]
	x = Matrix{ComplexF64}([0 im; -im 0])
	hop = Ω .* x
	z = [-1 0; 0 1]
	Is = one(x)
	Ib = one(zeros(d, d))

	jumpops = [randn(ComplexF64, 2, 2), randn(ComplexF64, 2, 2)]
	model = ImpurityLindbladian(lindbladoperator(hop, jumpops))

	Hbarebath = bosondensityoperator(d=d)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib) + kron(Is, Hbarebath) + kron(z, a' + a)

	jumpops2 = [kron(jump, Ib) for jump in jumpops]
	Lop = lindbladoperator(H, jumpops2)

	bs = NonAdditiveHyb(z)
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	# mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI = hybriddynamics_naive(lattice, corr, bs, trunc=trunc)
	# @test distance(mpsI, mpsI′) / norm(mpsI′) < tol
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mps = mult!(mpsK, mpsI, trunc=trunc)

	ρimp = _rand_dm(2)


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
	ids2N = [k for k in sampleidx(N) if k > 1]
	idsallN = [1; ids2N]
	for i in ids2N
		ind2 = ContourIndex(i, branch=:+)
		m = ContourOperator([ind2,ind1], [op, op])
		mps2 = apply!(m, lattice, deepcopy(mps))
		mps2 = initialstate!(mps2, lattice, ρimp)
		v = integrate(lattice, mps2) / Zval
		push!(corrs, v)
	end
	
	A = kron(op, Ib)
	corrs2 = correlation_2op_1t(Lop, A, A, ρ, 0:δt:t, reverse = false)
	corrs2 = corrs2[idsallN]
	@test norm(corrs - corrs2) / norm(corrs2) < tol


	# off-diagonal observables

	op1 = [0 0.8; 0 0]
	op2 = [0 0; 0.7*im 0]

	A1 = kron(op1, Ib)
	A2 = kron(op2, Ib)


	c1 = ContourIndex(1, branch=:+)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = apply!(ct, lattice, deepcopy(mps))
	mps2 = initialstate!(mps2, lattice, ρimp)
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	c2 = ContourIndex(1, branch=:+)
	for i in ids2N
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = apply!(ct, lattice, deepcopy(mps))
		mps2 = initialstate!(mps2, lattice, ρimp)
		v = integrate(lattice, mps2) / Zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, A1, A2, ρ, 0:δt:t, reverse = false)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

	op1 = [0 0.8*im; 0 0]
	op2 = [0 0; 0.7 0]


	A1 = kron(op1, Ib)
	A2 = kron(op2, Ib)

	c1 = ContourIndex(1, branch=:-)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = apply!(ct, lattice, deepcopy(mps))
	mps2 = initialstate!(mps2, lattice, ρimp)
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	for i in ids2N
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = apply!(ct, lattice, deepcopy(mps))
		mps2 = initialstate!(mps2, lattice, ρimp)
		v = integrate(lattice, mps2) / Zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, A1, A2, ρ, 0:δt:t, reverse = true)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


end
println("------------------------------------")
println("|     Quantum channel and rdm      |")
println("------------------------------------")

@testset "Rabi model: quantum channel and rdm" begin
	Ω = 0.5
	N = 10
	δt = 0.05
	t = N * δt
	β = 2
	chi = 100
	d = 50
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(N = N, δt=δt, contour=:real)

	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]
	hop = Ω .* z
	Is = one(x)
	Ib = one(zeros(d, d))
	model = ImpurityHamiltonian(hop)

	Hbarebath = bosondensityoperator(d=d)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib) + kron(Is, Hbarebath) + kron(y, a' + a)

	bs = NonAdditiveHyb(y)
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	mpsI = hybriddynamics_naive(lattice, corr, bs, trunc=trunc)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mps = mult!(mpsK, mpsI, trunc=trunc)

	# quantum channel against the exact full-system evolution
	Uop = exp((-im .* t) .* H)
	ρbath = exp(-β .* Hbarebath)
	Uop4 = reshape(Uop, (d, 2, d, 2))
	Uop4t = reshape(Uop', (d, 2, d, 2))

	@tensor map1[2, 6, 4, 5] := Uop4[1, 2, 3, 4] * ρbath[3, 7] * Uop4t[7, 5, 1, 6]
	map1 = reshape(map1, 4, 4)
	map2 = quantummap(lattice, mps)
	map2 = reshape(map2, 4, 4)
	map1 ./= tr(map1)
	map2 ./= tr(map2)
	@test distance(map1, map2) / norm(map1) < 2 * tol

	# reduced density matrix against the exact partial trace
	ρ1 = zeros(2, 2)
	ρ1[1, 1] = 1
	for ρimp in [ρ1, _rand_dm(2)]
		tmp = initialstate!(deepcopy(mps), lattice, ρimp)
		Zval = integrate(lattice, tmp)

		ρ = kron(ρimp, ρbath)
		ρout = Uop * ρ * Uop'
		ρout4 = reshape(ρout, (d, 2, d, 2))
		@tensor ρimpout[2, 3] := ρout4[1, 2, 1, 3]
		ρimpout ./= tr(ρimpout)
		ρimpout2 = rdm(lattice, tmp)
		ρimpout2 ./= tr(ρimpout2)
		@test distance(ρimpout, ρimpout2) / norm(ρimpout) < tol
	end
end
