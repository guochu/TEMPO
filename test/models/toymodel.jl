println("------------------------------------")
println("|             Toy Model            |")
println("------------------------------------")

# H = Ω*σy	+ σz(a + a†) + σz(b + b†) + a† a + b† b
@testset "Toy model: imaginary-time" begin

	Ω = 0.5
	N = 20
	δτ = 0.1
	β = N * δτ
	chi = 100
	d = 20
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(N = N, δτ=δτ, contour=:imag)

	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]
	hop = Ω .* y
	Is = one(x)
	Ib = one(zeros(d, d))
	model = BosonicImpurity(hop)

	nop = bosondensityoperator(d=d)
	Hbarebath = kron(nop, Ib) + kron(Ib, nop)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib, Ib) + kron(Is, Hbarebath) + kron(z, a' + a, Ib) + kron(z, Ib, a' + a)


	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	
	bs = AdditiveHyb([z[i,i] for i in 1:size(z,1)])

	spec = DiracDelta(1)

	bath = bosonicbath(spec, β=β)

	corr = correlationfunction(bath, lattice)

	mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI = mult!(mpsI, deepcopy(mpsI), trunc=trunc)
	mps = mult!(mpsK, mpsI, trunc=trunc)


	ρ = exp(-β * H)
	z1 = integrate(mps)
	z2 = tr(ρ) / tr(exp(-β .* Hbarebath))
	@test abs(z1 - z2) / abs(z1) < tol


	## diagonal observables
	op = [-0.73 0; 0 0.5]
	zdiag = [op[i,i] for i in 1:size(z, 1)]

	pos1 = index(lattice, 1)
	t = ADTTerm(pos1, zdiag .* zdiag )
	mps2 = apply!(t, copy(mps))
	v = integrate(mps2) / z1

	corrs = [v]
	for i in 2:N
		pos2 = index(lattice, i)
		t = ADTTerm((pos2,pos1), (zdiag, zdiag))
		# t = ADTTerm((i,1), reshape(kron(zdiag, zdiag), 2, 2))
		mps2 = apply!(t, copy(mps))
		v = integrate(mps2) / z1
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
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / z1

	corrs = [v]
	c2 = ContourIndex(1)
	for i in 2:N
		c1 = ContourIndex(i)
		ct = ContourOperator([c1, c2], [op1, op2])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / z1
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
	δt = 0.025
	β = 2
	t = N * δt
	chi = 100
	d = 20
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(N = N, δt=δt, contour=:real)

	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]
	hop = Ω .* y
	Is = one(x)
	Ib = one(zeros(d, d))
	model = BosonicImpurity(hop)

	nop = bosondensityoperator(d=d)
	Hbarebath = kron(nop, Ib) + kron(Ib, nop)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib, Ib) + kron(Is, Hbarebath) + kron(z, a' + a, Ib) + kron(z, Ib, a' + a)


	bs = AdditiveHyb([z[i,i] for i in 1:size(z,1)])
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI = mult!(mpsI, deepcopy(mpsI), trunc=trunc)

	ρ1 = zeros(2,2)
	ρ1[1,1] = 1
	ρ2 = 0.5 .* one(ρ1)

	for ρimp in [ρ1, ρ2, _rand_dm(2)]

		mpsK = sysdynamics(lattice, model, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
		mps = mult!(mpsK, mpsI, trunc=trunc)

		Zval = integrate(mps)

		
		ρ = kron(ρimp, exp(-β * Hbarebath)) 

		## diagonal observables
		op = [-0.73 0; 0 0.5]
		zdiag = [op[i,i] for i in 1:size(z, 1)]

		pos1 = index(lattice, 1, branch=:+)
		m = ADTTerm(pos1, zdiag .* zdiag )
		mps2 = apply!(m, copy(mps))
		v = integrate(mps2) / Zval

		corrs = [v]
		for i in 2:N
			pos2 = index(lattice, i, branch=:+)
			m = ADTTerm((pos2,pos1), (zdiag, zdiag))
			mps2 = apply!(m, copy(mps))
			v = integrate(mps2) / Zval
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
		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / Zval

		corrs = [v]
		c2 = ContourIndex(1, branch=:+)
		for i in 2:N
			c1 = ContourIndex(i, branch=:+)
			ct = ContourOperator([c1, c2], [op1, op2])

			mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
			mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
			mps2 = mult!(mpsK, mpsI, trunc=trunc)
			v = integrate(mps2) / Zval

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
		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / Zval

		corrs = [v]
		for i in 2:N
			c2 = ContourIndex(i, branch=:+)
			ct = ContourOperator([c1, c2], [op1, op2])

			mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
			mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
			mps2 = mult!(mpsK, mpsI, trunc=trunc)
			v = integrate(mps2) / Zval

			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = true)
		corrs2 = corrs2[1:length(corrs)]

		@test norm(corrs - corrs2) / norm(corrs2) < tol

	end

end



@testset "Toy model: mixed-time" begin

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

	lattice = ADTLattice(Nt = Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed)

	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]
	hop = Ω .* y
	Is = one(x)
	Ib = one(zeros(d, d))
	model = BosonicImpurity(hop)

	nop = bosondensityoperator(d=d)
	Hbarebath = kron(nop, Ib) + kron(Ib, nop)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib, Ib) + kron(Is, Hbarebath) + kron(z, a' + a, Ib) + kron(z, Ib, a' + a)


	op1 = [0 0.8; 0 0]
	op2 = [0 0; 0.7 0]

	A1 = kron(op1, Ib, Ib)
	A2 = kron(op2, Ib, Ib)


	bs = AdditiveHyb([z[i,i] for i in 1:size(z,1)])
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI = mult!(mpsI, deepcopy(mpsI), trunc=trunc)

	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	mps = mult!(mpsK, mpsI, trunc=trunc)

	Zval = integrate(mps)


	ρ = exp(-β .* H)

	# off-diagonal observables

	c1 = ContourIndex(1, branch=:+)

	ct = ContourOperator(c1, op1 * op2)
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	c2 = ContourIndex(1, branch=:+)
	for i in 2:Nt
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / Zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = false)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	c1 = ContourIndex(1, branch=:-)

	ct = ContourOperator(c1, op2 * op1)
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / Zval

	corrs = [v]
	for i in 2:Nt
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op2, op1])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / Zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(H, A2, A1, ρ, 0:δt:t, reverse = true)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end
