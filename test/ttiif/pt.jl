println("------------------------------------")
println("|           TTI-IF PT              |")
println("------------------------------------")


@testset "InfluenceOperator-NonAdditiveHyb: imaginary-time" begin
	β = 1
	δτ = 0.2
	tol = 1.0e-6
	N = round(Int, β/δτ)
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=50, ϵ=1.0e-6, add_back=0)
	algexpan = PronyExpansion(n=20, tol=1.0e-8)
	d = 2

	lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)
	corr = correlationfunction(bath, lattice)
	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]


	hyb = NonAdditiveHyb(y)
	op = hyb.op
	z2 = op * op
	zz = kron(op, op)

	mpo1 = influenceoperator(lattice, corr, hyb, algexpan=algexpan)

	local mpo2
	for i in 1:lattice.N, j in 1:lattice.N
		ind1 = ContourIndex(i)
		ind2 = ContourIndex(j)
		coef = index(corr, i, j)
		if coef != 0
			if ind1 == ind2
				m = coef .* z2
				t = FockTermS(lattice[ind1], m)
			else
				m = coef .* zz
				t = FockTermS((lattice[ind1], lattice[ind2]), reshape(m, (d,d,d,d)))
			end
			if @isdefined mpo2
				mpo2 += apply!(t, vacuumstate(lattice)) 
			else
				mpo2 = apply!(t, vacuumstate(lattice)) 
			end

		end
	end
	canonicalize!(mpo2, alg=Orthogonalize(trunc=trunc))
	@test distance(mpo1, mpo2) / norm(mpo2) < tol

	dt = 0.01

	mpo1 = dt * mpo1 + vacuumstate(lattice)
	mpo2 = influenceoperatorexponential(lattice, corr, dt, hyb, WII(), algexpan=algexpan)

	@test distance(mpo1, mpo2) / norm(mpo2) < dt

	# for algmult in (SVDCompression(D=50), DMRGMult1(trunc=truncdimcutoff(D=50,ϵ=1.0e-6)), DMRGMult2(trunc=truncdimcutoff(D=50,ϵ=1.0e-6)))
	# 	mps1 = differentialinfluencefunctional(lattice, corr, dt, WII(), algmult, band=band, algexpan=algexpan)
	# 	_n = norm(mps1)
	# 	mps2 = differentialinfluencefunctional(lattice, corr, dt, WI(), algmult, band=band, algexpan=algexpan)
	# 	mps3 = differentialinfluencefunctional(lattice, corr, dt, ComplexStepper(WI()), algmult, band=band, algexpan=algexpan)
	# 	mps4 = differentialinfluencefunctional(lattice, corr, dt, ComplexStepper(WII()), algmult, band=band, algexpan=algexpan)
	# 	@test distance(mps1, mps0) / _n < dt
	# 	@test distance(mps1, mps0) / _n < dt
	# 	@test distance(mps1, mps0) / _n < dt
	# 	@test distance(mps1, mps0) / _n < dt
	# end

end
