println("------------------------------------")
println("|    Independentbosons Model       |")
println("------------------------------------")



@testset "Independentbosons model: imaginary-time" begin
	δτ=0.1
	N = 10
	β = N * δτ
	chi = 50
	tol = 1.0e-2

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	trunc2 = truncdimcutoff(D=2*chi, ϵ=1.0e-10, add_back=0)

	lattice = PTLattice(N=N, δτ=δτ, contour=:imag)	

	for ϵ_d in (-1, 0, 1)

		nop = [1 0; 0 0]
		hyb = NonAdditiveHyb(nop)

		spec = Leggett(d=3, ωc=1)

		bath = bosonicbath(spec, β=β)
		corr = correlationfunction(bath, lattice)

		# mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)
		mpsI = hybriddynamics_naive(lattice, corr, hyb, trunc=trunc)
		# @test distance(mpsI, mpsI′) / norm(mpsI′) < tol


		model = BosonicImpurity(ϵ_d .* nop)
		mpsK = sysdynamics(lattice, model, trunc=trunc)

		adt = mult(mpsK, mpsI, trunc=trunc2)

		Zval = integrate(lattice, adt)

		# Green's function
		op1 = [0 0; 1 0]
		op2 = [0 1;0 0 ]

		c1 = ContourIndex(1)

		ct = ContourOperator(c1, op1 * op2)
		adt2 = apply!(ct, lattice, deepcopy(adt))
		v = integrate(lattice, adt2) / Zval

		corrs = [v]
		c2 = ContourIndex(1)
		for i in 2:N
			c1 = ContourIndex(i)
			ct = ContourOperator([c1, c2], [op1, op2])

			adt2 = apply!(ct, lattice, deepcopy(adt))
			v = integrate(lattice, adt2) / Zval
			push!(corrs, v)
		end

		corrs2 = independentbosons_Gτ(spec, β=β, ϵ_d=ϵ_d, Nτ=N)
		corrs2 = corrs2[1:length(corrs)]

		@test norm(corrs - corrs2) / norm(corrs2) < tol		

	end

end


@testset "Independentbosons model: mixed-time" begin
	δτ=0.1
	Nτ = 10
	β = Nτ * δτ
	δt = 0.02
	Nt = 10
	t = Nt * δt
	chi = 50
	tol = 1.0e-2

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	trunc2 = truncdimcutoff(D=2*chi, ϵ=1.0e-10, add_back=0)

	lattice = PTLattice(Nτ=Nτ, δτ=δτ, Nt=Nt, δt=δt, contour=:mixed)	

	ϵ_d = 0.5

	nop = [1 0; 0 0]
	hyb = NonAdditiveHyb(nop)

	spec = Leggett(d=3, ωc=1)

	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)

	# mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)
	mpsI = hybriddynamics_naive(lattice, corr, hyb, trunc=trunc)
	# @test distance(mpsI, mpsI′) / norm(mpsI′) < tol


	model = BosonicImpurity(ϵ_d .* nop)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	adt = mult(mpsK, mpsI, trunc=trunc2)

	Zval = integrate(lattice, adt)

	# greater Green's function
	op1 = [0 0; 1 0]
	op2 = [0 1;0 0 ]


	c1 = ContourIndex(1, branch=:+)

	ct = ContourOperator(c1, op1 * op2)
	mps2 = apply!(ct, lattice, deepcopy(adt))
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	c2 = ContourIndex(1, branch=:+)
	for i in 2:Nt
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mps2 = apply!(ct, lattice, deepcopy(adt))
		v = integrate(lattice, mps2) / Zval

		push!(corrs, v)
	end
	corrs = -im .* corrs

	corrs2 = [independentbosons_greater(spec, tj, β=β, ϵ_d=ϵ_d) for tj in 0:δt:t]
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

	# lesser Green's function
	c1 = ContourIndex(1, branch=:-)

	ct = ContourOperator(c1, op2 * op1)
	mps2 = apply!(ct, lattice, deepcopy(adt))
	v = integrate(lattice, mps2) / Zval

	corrs = [v]
	for i in 2:Nt
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op2, op1])

		mps2 = apply!(ct, lattice, deepcopy(adt))
		v = integrate(lattice, mps2) / Zval

		push!(corrs, v)
	end

	corrs = im .* corrs

	corrs2 = [independentbosons_lesser(spec, tj, β=β, ϵ_d=ϵ_d) for tj in 0:δt:t]
	corrs2 = corrs2[1:length(corrs)]

	# println(corrs)
	# println(corrs2)

	@test norm(corrs - corrs2) / norm(corrs2) < tol
end