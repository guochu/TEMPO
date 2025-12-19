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
	model = BosonicImpurity(xop)
	mps = sysdynamics(lattice, model, trunc=trunc)
	mps = boundarycondition!(mps, lattice, trunc=trunc)

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
	for i in 2:N
		pos2 = index(lattice, i)
		t = ADTTerm((pos2,pos1), (zdiag, zdiag))
		# t = ADTTerm((i,1), reshape(kron(zdiag, zdiag), 2, 2))
		mps2 = apply!(t, copy(mps))
		v = integrate(mps2) / integrate(mps)
		push!(corrs, v)
	end

	corrs2 = correlation_2op_1τ(xop, z, z, 0:δτ:β, β=β)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

	# off-diagonal observables
	op1 = [0 0; 1 0]
	op2 = [0 1; 0 0]

	c1 = ContourIndex(1)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
	mps2 = boundarycondition!(mps2, lattice, trunc=trunc)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	
	for i in 2:N
		c2 = ContourIndex(i)
		ct = ContourOperator([c2, c1], [op2, op1])

		mps2 = sysdynamics(lattice, model, ct, trunc=trunc)
		mps2 = boundarycondition!(mps2, lattice, trunc=trunc)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1τ(xop, op1, op2, 0:δτ:β, β=β)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end