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

	lattice = ADTLattice(N=N, δτ=δτ, contour=:imag)	

	ϵ_d = 0.5

	nop = [1 0; 0 0]
	hyb = AdditiveHyb([1, 0])

	spec = Leggett(d=3, ωc=1)

	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)

	mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)
	mpsI′ = hybriddynamics_naive(lattice, corr, hyb, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI′) < tol


	model = BosonicImpurity(nop)
	mpsK = sysdynamics(lattice, model, trunc=trunc)

	mpsK = boundarycondition!(mpsK, lattice, trunc=trunc)

	adt = mult(mpsK, mpsI, trunc=trunc2)

	Zval = integrate(adt)

	# Green's function
	op1 = [0 0; 1 0]
	op2 = [0 1;0 0 ]

	c1 = ContourIndex(1)

	ct = ContourOperator(c1, op1 * op2)
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, trunc=trunc)
	adt2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(adt2) / Zval

	corrs = [v]
	for i in 2:N
		c2 = ContourIndex(i)
		ct = ContourOperator([c2, c1], [op2, op1])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, trunc=trunc)
		adt2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(adt2) / Zval
		push!(corrs, v)
	end

	corrs2 = independentbosons_Gτ(spec, β=β, ϵ_d=1, Nτ=N)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end