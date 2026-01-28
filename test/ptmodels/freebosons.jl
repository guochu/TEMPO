println("------------------------------------")
println("|            Free Bosons           |")
println("------------------------------------")

omic_spectrum(w, α, wc) = 2α * w * exp(-w/wc)

@testset "Freebosons: imaginary-time" begin
	δτ=0.1
	N = 10
	β = N * δτ
	d = 10
	dw = 0.01
	α = 0.1
	wc = 2.
	chi = 30
	tol = 2.0e-2

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)	

	ϵ_d = 1.


	spec = spectrum(w -> omic_spectrum(w, α, wc), lb=0, ub=wc)

	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)


	adag = bosonadagoperator(d=d)
	hyb = NonDiagonalHyb(adag)
	algmult = DMRGMult1(trunc)
	algexpan = PronyExpansion(n=20, tol=1.0e-8)
	alg = TranslationInvariantIF(k=5, fast=true, algmult=algmult, algexpan=algexpan)

	mpsI = hybriddynamics(lattice, corr, hyb, alg)

	nop = bosondensityoperator(d=d)
	model = ImpurityHamiltonian(ϵ_d .* nop)
	mpsK = sysdynamics(lattice, model, trunc=trunc)

	adt = mpsK * mpsI

	Zval = integrate(lattice, adt)

	# Green's function
	op1 = adag'
	op2 = adag

	c1 = ContourIndex(1)

	ct = ContourOperator(c1, op1 * op2)
	adt2 = apply!(ct, lattice, copy(adt))
	v = integrate(lattice, adt2) / Zval

	corrs = [v]
	c2 = ContourIndex(1)
	for i in 2:N
		c1 = ContourIndex(i)
		ct = ContourOperator([c1, c2], [op1, op2])

		adt2 = apply!(ct, lattice, copy(adt))
		v = integrate(lattice, adt2) / Zval
		push!(corrs, v)
	end

	b2 = discretebath(bath, δw=dw)
	exact_model = Toulouse(b2, ϵ_d=ϵ_d)

	corrs2 = toulouse_Gτ(exact_model, 0:δτ:β)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol		
end


@testset "Freebosons: real-time" begin
	δt=0.1
	N = 10
	t = δt * N
	β = 2.
	d = 4
	dw = 0.01
	α = 0.1
	wc = 2.
	chi = 30
	tol = 1.0e-2

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	lattice = PTLattice(N=N, δt=δt, d=d, contour=:real)	

	ϵ_d = 1.


	spec = spectrum(w -> omic_spectrum(w, α, wc), lb=0, ub=wc)

	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)


	adag = bosonadagoperator(d=d)
	hyb = NonDiagonalHyb(adag)
	algmult = DMRGMult1(trunc)
	algexpan = PronyExpansion(n=20, tol=1.0e-8)
	alg = TranslationInvariantIF(k=5, fast=true, algmult=algmult, algexpan=algexpan)

	mpsI = hybriddynamics(lattice, corr, hyb, alg)

	nop = bosondensityoperator(d=d)
	model = ImpurityHamiltonian(ϵ_d .* nop)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	ρimp = bosonoccupationoperator(1, d=d)

	adt = mpsK * mpsI

	tmp = initialstate!(copy(adt), lattice, ρimp)
	Zval = integrate(lattice, adt)

	# Green's function
	op1 = adag'
	op2 = adag

	c1 = ContourIndex(1, :+)

	ct = ContourOperator(c1, op1 * op2)
	adt2 = apply!(ct, lattice, copy(adt))
	adt2 = initialstate!(adt2, lattice, ρimp)
	v = integrate(lattice, adt2) / Zval

	gt_corrs = [v]
	c2 = ContourIndex(1, :+)
	for i in 2:N
		c1 = ContourIndex(i, :+)
		ct = ContourOperator([c1, c2], [op1, op2])

		adt2 = apply!(ct, lattice, copy(adt))
		adt2 = initialstate!(adt2, lattice, ρimp)
		v = integrate(lattice, adt2) / Zval
		push!(gt_corrs, v)
	end
	gt_corrs = -im * gt_corrs

	c1 = ContourIndex(1, :-)
	ct = ContourOperator(c1, op2 * op1)
	adt2 = apply!(ct, lattice, copy(adt))
	adt2 = initialstate!(adt2, lattice, ρimp)
	v = integrate(lattice, adt2) / Zval

	lt_corrs = [v]	
	for i in 2:N
		c2 = ContourIndex(i, :+)
		ct = ContourOperator([c1, c2], [op2, op1])

		adt2 = apply!(ct, lattice, copy(adt))
		adt2 = initialstate!(adt2, lattice, ρimp)
		v = integrate(lattice, adt2) / Zval
		push!(lt_corrs, v)
	end	
	lt_corrs = im * lt_corrs


	b2 = discretebath(bath, δw=dw)
	exact_model = Toulouse(b2, ϵ_d=ϵ_d)

	gt_corrs2, lt_corrs2 = toulouse_neq_greater_lesser(exact_model, 0:δt:t, nsys=1)
	gt_corrs2 = gt_corrs2[1:length(gt_corrs)]
	lt_corrs2 = lt_corrs2[1:length(lt_corrs)]

	@test norm(gt_corrs - gt_corrs2) / norm(gt_corrs2) < tol	
	@test norm(lt_corrs - lt_corrs2) / norm(lt_corrs2) < tol	
end