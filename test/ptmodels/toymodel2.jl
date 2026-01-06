println("------------------------------------")
println("|           Toy JC Model           |")
println("------------------------------------")

# H = Ω*σz	+ (A†a + Aa†) + 2a† a
@testset "Toy JC model: imaginary-time" begin

	Ω = 0.5
	N = 20
	δτ = 0.05
	β = N * δτ
	chi = 50
	d = 50
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(N = N, δτ=δτ, contour=:imag)

	p = spin_half_matrices()
	x, y, z, sp = p["x"], p["y"], p["z"], p["+"]
	hop = Ω .* z
	Is = one(x)
	Ib = one(zeros(d, d))
	model = BosonicImpurity(hop)

	Hbarebath = bosondensityoperator(d=d)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib) + kron(Is, Hbarebath) + kron(sp, a) + kron(sp', a')

	mpsK = sysdynamics(lattice, model, trunc=trunc)
	
	hyb = NonDiagonalHyb(sp)

	spec = DiracDelta(1)

	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)

	algmult = DMRGMult1(trunc)
	algexpan = PronyExpansion(n=20, tol=1.0e-8)
	alg = TranslationInvariantIF(k=5, fast=true, algmult=algmult, algexpan=algexpan)
	mpsI = hybriddynamics(lattice, corr, hyb, alg)
	mps = mult!(mpsK, mpsI, trunc=trunc)

	Zval = integrate(lattice, mps)

	ρ = exp(-β * H)
	Zval2 = tr(ρ) / tr(exp(-β .* Hbarebath))
	# println(Zval, " ", Zval2)
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
	

	A = kron(op, Ib)

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
	

	A1 = kron(op1, Ib)
	A2 = kron(op2, Ib)

	corrs2 = correlation_2op_1τ(H, A1, A2, 0:δτ:β, β=β)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end

@testset "Toy JC model: real-time" begin

	Ω = 0.5
	N = 10
	δt = 0.05
	β = 2
	t = N * δt
	chi = 50
	d = 50
	tol = 2.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(N = N, δt=δt, contour=:real)

	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]
	hop = Ω .* z
	Is = one(x)
	Ib = one(zeros(d, d))
	model = BosonicImpurity(hop)
	sp = randn(ComplexF64, 2, 2)
	sp ./= norm(sp)

	Hbarebath = bosondensityoperator(d=d)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib) + kron(Is, Hbarebath) + kron(sp, a) + kron(sp', a')

	bs = NonDiagonalHyb(sp)
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)

	algmult = DMRGMult1(trunc)
	algexpan = PronyExpansion(n=20, tol=1.0e-8)
	alg = TranslationInvariantIF(k=5, fast=true, algmult=algmult, algexpan=algexpan, verbosity=2)
	mpsI = hybriddynamics(lattice, corr, bs, alg)
	# @test distance(mpsI, mpsI′) / norm(mpsI′) < tol
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mps = mult!(mpsK, mpsI, trunc=trunc)

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
		
		A = kron(op, Ib)
		corrs2 = correlation_2op_1t(H, A, A, ρ, 0:δt:t, reverse = false)
		corrs2 = corrs2[1:length(corrs)]
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


		A1 = kron(op1, Ib)
		A2 = kron(op2, Ib)

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
