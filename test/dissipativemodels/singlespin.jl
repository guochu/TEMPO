println("------------------------------------")
println("|      Dissipative Single Spin     |")
println("------------------------------------")


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
	for i in 2:N
		pos2 = index(lattice, i, branch=:+)
		m = ADTTerm((pos2,pos1), (zdiag, zdiag))
		mps2 = apply!(m, deepcopy(mps))
		v = integrate(mps2) / integrate(mps)
		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, z, z, ρimp, 0:δt:t, reverse = false)

	corrs2 = corrs2[1:length(corrs)]

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
	for i in 2:N
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
		mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, op1, op2, ρimp, 0:δt:t, reverse = false)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	c1 = ContourIndex(1, branch=:-)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
	mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	for i in 2:N
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
		mps2 = boundarycondition!(mps2, lattice, ρ₀=ρimp)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, op1, op2, ρimp, 0:δt:t, reverse = true)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


end
