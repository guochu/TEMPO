println("------------------------------------")
println("|           Single Spin            |")
println("------------------------------------")


@testset "Single spin: imaginary-time" begin
	Ω = 0.5
	N = 10
	δτ = 0.1
	β = N * δτ
	chi = 50
	tol = 1.0e-6
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(N = N, δτ=δτ, contour=:imag)

	xop = Ω .* [0 1; 1 0]
	model = ImpurityHamiltonian(xop)
	mps = sysdynamics(lattice, model, trunc=trunc)
	mps = boundarycondition!(mps, lattice)

	ρ = exp(-β * xop)
	@test tr(ρ) ≈ integrate(mps)

	z = [1 0; 0 0]
	zdiag = [z[i,i] for i in 1:size(z, 1)]
	

	pos1 = index(lattice, 1)
	t = ADTTerm((pos1,), zdiag .* zdiag )

	mps2 = apply!(t, copy(mps))

	v = integrate(mps2) / integrate(mps)
	@test v ≈ tr(z * ρ) / tr(ρ)

	corrs = [v]
	ids2N = [k for k in sampleidx(N) if k > 1]
	idsallN = [1; ids2N]
	for i in ids2N
		pos2 = index(lattice, i)
		t = ADTTerm((pos2,pos1), (zdiag, zdiag))
		# t = ADTTerm((i,1), reshape(kron(zdiag, zdiag), 2, 2))
		mps2 = apply!(t, copy(mps))
		v = integrate(mps2) / integrate(mps)
		push!(corrs, v)
	end

	corrs2 = correlation_2op_1τ(xop, z, z, 0:δτ:β, β=β)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

	# off-diagonal observables
	op1 = [0 0; 1 0]
	op2 = [0 1; 0 0]

	c1 = ContourIndex(1)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
	mps2 = boundarycondition!(mps2, lattice)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	
	for i in ids2N
		c2 = ContourIndex(i)
		ct = ContourOperator([c2, c1], [op2, op1])

		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
		mps2 = boundarycondition!(mps2, lattice)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1τ(xop, op1, op2, 0:δτ:β, β=β)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end


@testset "Single spin: real-time" begin
	Ω = 0.5
	N = 10
	δt = 0.1
	t = N * δt
	β = 1
	chi = 50
	tol = 1.0e-6
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(N = N, δt=δt, contour=:real)

	ρ1 = zeros(2,2)
	ρ1[1,1] = 1
	ρ2 = 0.5 .* one(ρ1)

	z = [0.6 0; 0 0.3]
	zdiag = [z[i,i] for i in 1:size(z, 1)]
	
	xop = Ω .* [0 1; 1 0]
	model = ImpurityHamiltonian(xop)

	for ρimp in [ρ1, ρ2]

		mps = sysdynamics(lattice, model, trunc=trunc)
		mps = boundarycondition!(mps, lattice, ρ₀=ρimp)

		pos1 = index(lattice, 1, branch=:+)
		m = ADTTerm((pos1,), zdiag .* zdiag )

		mps2 = apply!(m, deepcopy(mps))

		# diagonal observables
		v = integrate(mps2) / integrate(mps)
		@test v ≈ tr(z*z * ρimp) 

		corrs = [v]
		ids2N = [k for k in sampleidx(N) if k > 1]
		idsallN = [1; ids2N]
		for i in ids2N
			pos2 = index(lattice, i, branch=:+)
			m = ADTTerm((pos2,pos1), (zdiag, zdiag))
			mps2 = apply!(m, deepcopy(mps))
			v = integrate(mps2) / integrate(mps)
			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(xop, z, z, ρimp, 0:δt:t, reverse = false)

		corrs2 = corrs2[idsallN]

		@test norm(corrs - corrs2) / norm(corrs2) < tol

		# off-diagonal observables
		op1 = [0 0.8; 0 0]
		op2 = [0 0; 0.7 0]

		c1 = ContourIndex(1, branch=:+)

		ct = ContourOperator(c1, op1 * op2)
		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
		mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
		v = integrate(mps2) / integrate(mps)

		corrs = [v]
		c2 = ContourIndex(1, branch=:+)
		for i in ids2N
			c1 = ContourIndex(i, branch=:+)
			ct = ContourOperator([c1, c2], [op1, op2])

			mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
			mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
			v = integrate(mps2) / integrate(mps)

			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(xop, op1, op2, ρimp, 0:δt:t, reverse = false)
		corrs2 = corrs2[idsallN]

		@test norm(corrs - corrs2) / norm(corrs2) < tol


		c1 = ContourIndex(1, branch=:-)

		ct = ContourOperator(c1, op1 * op2)
		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
		mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
		v = integrate(mps2) / integrate(mps)

		corrs = [v]
		for i in ids2N
			c2 = ContourIndex(i, branch=:+)
			ct = ContourOperator([c1, c2], [op1, op2])

			mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
			mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
			v = integrate(mps2) / integrate(mps)

			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(xop, op1, op2, ρimp, 0:δt:t, reverse = true)
		corrs2 = corrs2[idsallN]

		@test norm(corrs - corrs2) / norm(corrs2) < tol

	end

end


@testset "Single spin: mixed-time" begin
	Ω = 0.5
	Nt = 5
	δt = 0.05
	t = Nt * δt
	Nτ = 10
	δτ = 0.1
	β = Nτ * δτ
	chi = 100
	tol = 1.0e-6
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(Nτ = Nτ, δτ=δτ, Nt=Nt, δt=δt, contour=:mixed)

	
	xop = Ω .* [0 1; 1 0]
	model = ImpurityHamiltonian(xop)
	ρ = exp(-β .* xop)


	mps = sysdynamics(lattice, model, trunc=trunc)
	mps = boundarycondition!(mps, lattice)


	# off-diagonal observables
	op1 = [0 0.8; 0 0]
	op2 = [0 0; 0.7 0]

	c1 = ContourIndex(1, branch=:+)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
	mps2 = boundarycondition!(mps2, lattice)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	c2 = ContourIndex(1, branch=:+)
	ids2Nt = [k for k in sampleidx(Nt) if k > 1]
	idsallNt = [1; ids2Nt]
	for i in ids2Nt
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
		mps2 = boundarycondition!(mps2, lattice)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(xop, op1, op2, ρ, 0:δt:t, reverse = false)
	corrs2 = corrs2[idsallNt]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	c1 = ContourIndex(1, branch=:-)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
	mps2 = boundarycondition!(mps2, lattice)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	for i in ids2Nt
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
		mps2 = boundarycondition!(mps2, lattice)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(xop, op1, op2, ρ, 0:δt:t, reverse = true)
	corrs2 = corrs2[idsallNt]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


end
println("------------------------------------")
println("|      Dissipative Single Spin     |")
println("------------------------------------")


@testset "Single spin: real-time (dissipative)" begin
	Ω = 0.5
	N = 10
	δt = 0.1
	t = N * δt
	β = 1
	chi = 50
	tol = 1.0e-6
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(N = N, δt=δt, contour=:real)

	d = 2
	ρimp = _rand_dm(d)
	# ρimp = [0.7 0; 0 0.3]

	Lop = _rand_lindblad(d)
	# Lop = lindbladoperator(_rand_ham(ComplexF64, d), [])

	model = ImpurityLindbladian(Lop)

	z = [1 0; 0 0]
	zdiag = [z[i,i] for i in 1:size(z, 1)]


	mps = sysdynamics(lattice, model, trunc=trunc)
	mps = boundarycondition!(mps, lattice, ρ₀=ρimp)

	pos1 = index(lattice, 1, branch=:+)
	m = ADTTerm((pos1,), zdiag .* zdiag )

	mps2 = apply!(m, deepcopy(mps))

	# diagonal observables
	v = integrate(mps2) / integrate(mps)
	@test v ≈ tr(z*z * ρimp) 

	corrs = [v]
	ids2N = [k for k in sampleidx(N) if k > 1]
	idsallN = [1; ids2N]
	for i in ids2N
		pos2 = index(lattice, i, branch=:+)
		m = ADTTerm((pos2,pos1), (zdiag, zdiag))
		mps2 = apply!(m, deepcopy(mps))
		v = integrate(mps2) / integrate(mps)
		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, z, z, ρimp, 0:δt:t, reverse = false)

	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

	# off-diagonal observables
	op1 = [0 0.8; 0 0]
	op2 = [0 0; 0.7 0]

	c1 = ContourIndex(1, branch=:+)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
	mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	c2 = ContourIndex(1, branch=:+)
	for i in ids2N
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
		mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, op1, op2, ρimp, 0:δt:t, reverse = false)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	c1 = ContourIndex(1, branch=:-)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
	mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	for i in ids2N
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
		mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, op1, op2, ρimp, 0:δt:t, reverse = true)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


end
