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

	lattice = PTLattice(N = N, δτ=δτ, contour=:imag)

	xop = Ω .* [0 1; 1 0]
	model = BosonicImpurity(xop)
	mps = sysdynamics(lattice, model, trunc=trunc)

	ρ = exp(-β * xop)
	@test tr(ρ) ≈ integrate(lattice, mps)

	z = [1 0; 0 0]
	

	ind1 = ContourIndex(1)
	t = ContourOperator(ind1, z  )

	mps2 = apply!(t, lattice, copy(mps))

	v = integrate(lattice, mps2) / integrate(lattice, mps)
	@test v ≈ tr(z * ρ) / tr(ρ)


	# off-diagonal observables
	op1 = [0 0; 1 0]
	op2 = [0 1; 0 0]

	ind1 = ContourIndex(1)

	t = ContourOperator(ind1, op1 * op2)
	mps2 = apply!(t, lattice, copy(mps))
	v = integrate(lattice, mps2) / integrate(lattice, mps)

	corrs = [v]
	
	for i in 2:N
		ind2 = ContourIndex(i)

		t = ContourOperator([ind2, ind1], [op2, op1])
		mps2 = apply!(t, lattice, copy(mps))
		v = integrate(lattice, mps2) / integrate(lattice, mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1τ(xop, op1, op2, 0:δτ:β, β=β)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end


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

	ρ1 = zeros(2,2)
	ρ1[1,1] = 1
	ρ2 = 0.5 .* one(ρ1)

	z = [0.6 0; 0 0.3]
	
	# hop = Ω .* [0 1; 1 0]
	hop = Ω .* pauli_y()
	model = BosonicImpurity(hop)

	for ρimp in [ρ1,]

		mps = sysdynamics(lattice, model, trunc=trunc)
		# mps = initialstate!(mps, lattice, ρimp, trunc=trunc)

		ind1 = ContourIndex(1, branch=:+)
		m = ContourOperator(ind1, z * z )

		mps2 = apply!(m, lattice, deepcopy(mps), aheads=true)
		mps2 = initialstate!(mps2, lattice, ρimp, trunc=trunc)
		tmp = initialstate!(deepcopy(mps), lattice, ρimp, trunc=trunc)

		Zval = integrate(lattice, tmp)
		v = integrate(lattice, mps2) / Zval
		@test v ≈ tr(z*z * ρimp) 


		# off-diagonal observables
		op1 = [0 0.8; 0 0]
		op2 = [0 0; 0.7 0]

		ind1 = ContourIndex(1, branch=:+)

		m = ContourOperator(ind1, op1 * op2)
		mps2 = apply!(m, lattice, deepcopy(mps), aheads=true)
		mps2 = initialstate!(mps2, lattice, ρimp, trunc=trunc)
		v = integrate(lattice, mps2) / Zval

		corrs = [v]
		ind2 = ContourIndex(1, branch=:+)
		for i in 2:N
			ind1 = ContourIndex(i, branch=:+)
			m = ContourOperator([ind1, ind2], [op1, op2])

			mps2 = apply!(m, lattice, deepcopy(mps), aheads=true)
			mps2 = initialstate!(mps2, lattice, ρimp, trunc=trunc)
			v = integrate(lattice, mps2) / Zval

			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(hop, op1, op2, ρimp, 0:δt:t, reverse = false)
		corrs2 = corrs2[1:length(corrs)]

		@test norm(corrs - corrs2) / norm(corrs2) < tol


		# pos1 = index(lattice, 1, branch=:-) 
		# println("pos1 = ", pos1)

		# m = FockTerm(pos1, op1 * op2)
		# mps2 = apply!(m, deepcopy(mps), aheads=true)
		# mps2 = initialstate!(mps2, lattice, ρimp, trunc=trunc)
		# v = integrate(lattice, mps2) / Zval

		# corrs = [v]
		# for i in 2:N
		# 	pos2 = index(lattice, i, branch=:+)
		# 	println("pos2 = ", pos2)
		# 	# m = FockTerm([pos1, pos2], [Matrix(transpose(op1)), op2])
		# 	# mps2 = apply!(m, deepcopy(mps), aheads=(true, true))


		# 	m = ContourOperator([ContourIndex(1, :-), ContourIndex(i, :+)], [op1, op2])
		# 	mps2 = apply!(m, lattice, deepcopy(mps))
		# 	mps2 = initialstate!(mps2, lattice, ρimp, trunc=trunc)

		# 	v = integrate(lattice, mps2) / Zval

		# 	push!(corrs, v)
		# end

		# corrs2 = correlation_2op_1t(hop, op1, op2, ρimp, 0:δt:t, reverse = true)
		# corrs2 = corrs2[1:length(corrs)]

		# println(corrs)
		# println(corrs2)

		# @test norm(corrs - corrs2) / norm(corrs2) < tol

	end

end


# @testset "Single spin: mixed-time" begin
# 	Ω = 0.5
# 	Nt = 5
# 	δt = 0.05
# 	t = Nt * δt
# 	Nτ = 10
# 	δτ = 0.1
# 	β = Nτ * δτ
# 	chi = 100
# 	tol = 1.0e-6
# 	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

# 	lattice = ADTLattice(Nτ = Nτ, δτ=δτ, Nt=Nt, δt=δt, contour=:mixed)

	
# 	xop = Ω .* [0 1; 1 0]
# 	model = BosonicImpurity(xop)
# 	ρ = exp(-β .* xop)


# 	mps = sysdynamics(lattice, model, trunc=trunc)
# 	mps = boundarycondition!(mps, lattice, trunc=trunc)


# 	# off-diagonal observables
# 	op1 = [0 0.8; 0 0]
# 	op2 = [0 0; 0.7 0]

# 	c1 = ContourIndex(1, branch=:+)

# 	ct = ContourOperator(c1, op1 * op2)
# 	mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
# 	mps2 = boundarycondition!(mps2, lattice, trunc=trunc)
# 	v = integrate(mps2) / integrate(mps)

# 	corrs = [v]
# 	c2 = ContourIndex(1, branch=:+)
# 	for i in 2:Nt
# 		c1 = ContourIndex(i, branch=:+)
# 		ct = ContourOperator([c1, c2], [op1, op2])

# 		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
# 		mps2 = boundarycondition!(mps2, lattice, trunc=trunc)
# 		v = integrate(mps2) / integrate(mps)

# 		push!(corrs, v)
# 	end

# 	corrs2 = correlation_2op_1t(xop, op1, op2, ρ, 0:δt:t, reverse = false)
# 	corrs2 = corrs2[1:length(corrs)]

# 	@test norm(corrs - corrs2) / norm(corrs2) < tol


# 	c1 = ContourIndex(1, branch=:-)

# 	ct = ContourOperator(c1, op1 * op2)
# 	mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
# 	mps2 = boundarycondition!(mps2, lattice, trunc=trunc)
# 	v = integrate(mps2) / integrate(mps)

# 	corrs = [v]
# 	for i in 2:Nt
# 		c2 = ContourIndex(i, branch=:+)
# 		ct = ContourOperator([c1, c2], [op1, op2])

# 		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
# 		mps2 = boundarycondition!(mps2, lattice, trunc=trunc)
# 		v = integrate(mps2) / integrate(mps)

# 		push!(corrs, v)
# 	end

# 	corrs2 = correlation_2op_1t(xop, op1, op2, ρ, 0:δt:t, reverse = true)
# 	corrs2 = corrs2[1:length(corrs)]

# 	@test norm(corrs - corrs2) / norm(corrs2) < tol


# end