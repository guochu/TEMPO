println("------------------------------------")
println("|      Dissipative Single Spin     |")
println("------------------------------------")



@testset "Single spin: real-time" begin
	Ω = 0.5
	N = 10
	δt = 0.5
	t = N * δt
	β = 1
	chi = 50
	tol = 1.0e-6
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = PTLattice(N = N, δt=δt, contour=:real)

	d = 2
	ρimp = _rand_dm(d)
	# ρimp = [0.7 0; 0 0.3]

	Lop = _rand_lindblad(d)
	# Lop = lindbladoperator(_rand_ham(ComplexF64, d), [])
	
	# hop = Ω .* [0 1; 1 0]
	z = [1 0; 0 0.]
	model = DissipativeImpurity(Lop)

	mps = sysdynamics(lattice, model, trunc=trunc)
	# mps = initialstate!(mps, lattice, ρimp)

	ind1 = ContourIndex(1, branch=:+)
	m = ContourOperator(ind1, z * z )

	mps2 = apply!(m, lattice, deepcopy(mps), aheads=true)
	mps2 = initialstate!(mps2, lattice, ρimp)
	tmp = initialstate!(deepcopy(mps), lattice, ρimp)

	Zval = integrate(lattice, tmp)
	v = integrate(lattice, mps2) / Zval
	@test v ≈ tr(z*z * ρimp) 


	# off-diagonal observables
	op1 = [0 0.8; 0 0]
	op2 = [0 0; 0.7 0]

	ind1 = ContourIndex(1, branch=:+)

	m = ContourOperator(ind1, op1 * op2)
	mps2 = apply!(m, lattice, deepcopy(mps), aheads=true)
	mps2 = initialstate!(mps2, lattice, ρimp)
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	ind2 = ContourIndex(1, branch=:+)
	for i in 2:N
		ind1 = ContourIndex(i, branch=:+)
		m = ContourOperator([ind1, ind2], [op1, op2])

		mps2 = apply!(m, lattice, deepcopy(mps), aheads=true)
		mps2 = initialstate!(mps2, lattice, ρimp)
		v = integrate(lattice, mps2) / Zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, op1, op2, ρimp, 0:δt:t, reverse = false)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	ind1 = ContourIndex(1, branch=:-) 

	m = ContourOperator(ind1, op1 * op2)
	mps2 = apply!(m, lattice, deepcopy(mps), aheads=true)
	mps2 = initialstate!(mps2, lattice, ρimp)
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	for i in 2:N
		ind2 = ContourIndex(i, branch=:+)

		m = ContourOperator([ind1, ind2], [op1, op2])
		mps2 = apply!(m, lattice, deepcopy(mps))
		mps2 = initialstate!(mps2, lattice, ρimp)

		v = integrate(lattice, mps2) / Zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, op1, op2, ρimp, 0:δt:t, reverse = true)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end
