println("------------------------------------")
println("|     Dissipative Rabi Model       |")
println("------------------------------------")

@testset "Rabi model: real-time" begin

	Ω = 0.5
	N = 10
	δt = 0.05
	β = 2
	t = N * δt
	chi = 50
	d = 20
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(N = N, δt=δt, contour=:real)

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


	bs = AdditiveHyb([z[i,i] for i in 1:size(z,1)])
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI′ = hybriddynamics_naive(lattice, corr, bs, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI′) < tol


	ρimp = _rand_dm(2)



	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
	mps = mult!(mpsK, mpsI, trunc=trunc)

	
	ρ = kron(ρimp, exp(-β * Hbarebath)) 

	## diagonal observables
	op = [-0.73 0; 0 0.5]
	zdiag = [op[i,i] for i in 1:size(z, 1)]

	pos1 = index(lattice, 1, branch=:+)
	m = ADTTerm(pos1, zdiag .* zdiag )
	mps2 = apply!(m, copy(mps))
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	for i in 2:N
		pos2 = index(lattice, i, branch=:+)
		m = ADTTerm((pos2,pos1), (zdiag, zdiag))
		mps2 = apply!(m, copy(mps))
		v = integrate(mps2) / integrate(mps)
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
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	c2 = ContourIndex(1, branch=:+)
	for i in 2:N
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / integrate(mps)

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
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	for i in 2:N
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, A1, A2, ρ, 0:δt:t, reverse = true)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end

