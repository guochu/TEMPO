println("------------------------------------")
println("|     Dissipative Rabi Model       |")
println("------------------------------------")

@testset "Rabi model: real-time" begin

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

	# x = [0 1; 1 0]
	x = Matrix{ComplexF64}([0 im; -im 0])
	hop = Ω .* x
	z = [-1 0; 0 1]
	Is = one(x)
	Ib = one(zeros(d, d))

	jumpops = [randn(ComplexF64, 2, 2), randn(ComplexF64, 2, 2)]
	model = DissipativeImpurity(lindbladoperator(hop, jumpops))

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
	for i in 2:N
		ind2 = ContourIndex(i, branch=:+)
		m = ContourOperator([ind2,ind1], [op, op])
		mps2 = apply!(m, lattice, deepcopy(mps))
		mps2 = initialstate!(mps2, lattice, ρimp)
		v = integrate(lattice, mps2) / Zval
		push!(corrs, v)
	end
	
	A = kron(op, Ib)
	corrs2 = correlation_2op_1t(Lop, A, A, ρ, 0:δt:t, reverse = false)
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

	corrs2 = correlation_2op_1t(Lop, A1, A2, ρ, 0:δt:t, reverse = false)
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

	corrs2 = correlation_2op_1t(Lop, A1, A2, ρ, 0:δt:t, reverse = true)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


end