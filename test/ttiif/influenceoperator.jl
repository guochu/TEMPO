println("------------------------------------")
println("|    TTI-IF InfluenceOperator      |")
println("------------------------------------")


@testset "TTI-IF-InfluenceOperatot: imaginary-time" begin
	β = 1
	δτ = 0.2
	tol = 1.0e-6
	N = round(Int, β/δτ)
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=50, ϵ=1.0e-6, add_back=0)
	algexpan = PronyExpansion(n=20, tol=1.0e-8)
	d = 2

	lattice = ADTLattice(N=N, δτ=δτ, d=d, contour=:imag)
	corr = correlationfunction(bath, lattice)

	hyb = AdditiveHyb(randn(Float64, d))

	z = hyb.op
	z2 = z .* z
	zz = reshape(kron(z, z), d, d)

	mpo1 = influenceoperator(lattice, corr, hyb, algexpan=algexpan)

	orth = Orthogonalize(SVD(), trunc)
	local mpo2
	for i in 1:lattice.N, j in 1:lattice.N
		ind1, ind2 = ContourIndex(i+1), ContourIndex(j+1)
		coef = index(corr, i, j)
		if coef != 0
			if ind1 == ind2
				m = coef .* z2
				# m = coef .* op1 * op2 
				t = ADTTerm(lattice[ind1], m)
			else
				m = coef .* zz
				t = ADTTerm((lattice[ind1], lattice[ind2]), m)
			end

			# apply!(t, mpo2) 
			if @isdefined mpo2
				mpo2 += apply!(t, vacuumstate(lattice)) 
			else
				mpo2 = apply!(t, vacuumstate(lattice)) 
			end
			canonicalize!(mpo2, alg=orth)
		end
	end
	# canonicalize!(mpo2, alg=Orthogonalize(trunc=trunc))
	# println(distance(mpo1, mpo2), " ", norm(mpo1), " ", norm(mpo2))
	@test distance(mpo1, mpo2) / norm(mpo2) < tol

	dt = 0.01

	mpo1 = dt * mpo1 + vacuumstate(lattice)
	mps0 = influenceoperatorexponential(lattice, corr, dt, hyb, WII(), algexpan=algexpan)

	@test distance(mpo1, mps0) / norm(mps0) < dt

	for algmult in (SVDCompression(D=50), DMRGMult1(trunc=truncdimcutoff(D=50,ϵ=1.0e-6)))
		mps1 = differentialinfluencefunctional(lattice, corr, dt, hyb, WII(), algmult, algexpan=algexpan)
		_n = norm(mps1)
		mps2 = differentialinfluencefunctional(lattice, corr, dt, hyb, WI(), algmult, algexpan=algexpan)
		mps3 = differentialinfluencefunctional(lattice, corr, dt, hyb, ComplexStepper(WI()), algmult, algexpan=algexpan)
		mps4 = differentialinfluencefunctional(lattice, corr, dt, hyb, ComplexStepper(WII()), algmult, algexpan=algexpan)
		@test distance(mps1, mps0) / _n < dt
		@test distance(mps1, mps0) / _n < dt
		@test distance(mps1, mps0) / _n < dt
		@test distance(mps1, mps0) / _n < dt
	end

end

# @testset "TTI-IF-InfluenceOperatot: real-time" begin

# 	function __get_contour_op(lattice, ind1::ContourIndex, ind2::ContourIndex, z1::AbstractMatrix, z2::AbstractMatrix, coef)
# 		d = lattice.d
# 		if branch(ind1) == :-
# 			z1 = transpose(z1)
# 		end
# 		if branch(ind2) == :-
# 			z2 = transpose(z2)
# 		end
# 		pos1, pos2 = lattice[ind1], lattice[ind2]
# 		if pos1 == pos2
# 			m = (coef/2) .* (z1 * z2 + z2 * z1)
# 			t = FockTermS(pos1, m)
# 		else
# 			zz = kron(z2, z1)
# 			m = coef .* zz
# 			t = FockTermS((pos1, pos2), reshape(m, (d,d,d,d)))
# 		end
# 		return t
# 	end

# 	β = 1.
# 	δt = 0.1
# 	N = 8

# 	spec = Leggett(d=1, ωc=1)
# 	# spec = DiracDelta(1)
# 	bath = bosonicbath(spec, β=β, μ=0)
# 	trunc = truncdimcutoff(D=100, ϵ=1.0e-8, add_back=0)
# 	algexpan = PronyExpansion(n=20, tol=1.0e-8)
# 	d = 2
# 	tol = 1.0e-5



# 	lattice = PTLattice(N=N, δt=δt, contour=:real)
# 	corr = correlationfunction(bath, lattice)

# 	# hyb = NonAdditiveHyb(_rand_ham(ComplexF64, d))
# 	for hyb in (NonAdditiveHyb(_rand_ham(d)), NonDiagonalHyb(randn(ComplexF64, d, d)))

# 		op1, op2 = pairop(hyb)

# 		mps_pp, mps_pm, mps_mp, mps_mm = influenceoperator(lattice, corr, hyb, algexpan=algexpan)

# 		vc = vacuumstate(lattice)
# 		mps2_pp = nothing
# 		for i in 1:lattice.N, j in 1:lattice.N
# 			ind1, ind2 = ContourIndex(i, :+), ContourIndex(j, :+)

# 			coef = index(corr, i, j, b1=:+, b2=:+)
# 			t = __get_contour_op(lattice, ind1, ind2, op1, op2, coef)

# 			if !isnothing(mps2_pp)
# 				mps2_pp += apply!(t, vacuumstate(scalartype(hyb), lattice)) 
# 			else
# 				mps2_pp = apply!(t, vacuumstate(scalartype(hyb), lattice)) 
# 			end
# 		end

# 		# println(norm(mps_pp), " ", norm(mps2_pp), " ", dot(mps_pp, mps2_pp), " ", distance(mps_pp, mps2_pp))
# 		@test distance(mps_pp, mps2_pp) / norm(mps2_pp) < tol

# 		mps2_pm = nothing
# 		for i in 1:lattice.N, j in 1:lattice.N
# 			ind1, ind2 = ContourIndex(i, :+), ContourIndex(j, :-)
# 			coef = index(corr, i, j, b1=:+, b2=:-)
# 			t = __get_contour_op(lattice, ind1, ind2, op1, op2, coef)
# 			if !isnothing(mps2_pm)
# 				mps2_pm += apply!(t, vacuumstate(scalartype(hyb), lattice)) 
# 			else
# 				mps2_pm = apply!(t, vacuumstate(scalartype(hyb), lattice)) 
# 			end

# 		end
# 		# println(norm(mps_pm), " ", norm(mps2_pm), " ",  dot(mps_pm, mps2_pm), " ", distance(mps_pm, mps2_pm))
# 		@test distance(mps_pm, mps2_pm) / norm(mps2_pm) < tol

# 		mps2_mp = nothing
# 		for i in 1:lattice.N, j in 1:lattice.N
# 			ind1, ind2 = ContourIndex(i, :-), ContourIndex(j, :+)
# 			coef = index(corr, i, j, b1=:-, b2=:+)
# 			t = __get_contour_op(lattice, ind1, ind2, op1, op2, coef)
# 			if !isnothing(mps2_mp)
# 				mps2_mp += apply!(t, vacuumstate(scalartype(hyb), lattice)) 
# 			else
# 				mps2_mp = apply!(t, vacuumstate(scalartype(hyb), lattice)) 
# 			end

# 		end
# 		# println(norm(mps_mp), " ", norm(mps2_mp), " ", dot(mps_mp, mps2_mp), " ", distance(mps_mp, mps2_mp))
# 		@test distance(mps_mp, mps2_mp) / norm(mps2_mp) < tol				

# 		mps2_mm = nothing
# 		for i in 1:lattice.N, j in 1:lattice.N
# 			ind1, ind2 = ContourIndex(i, :-), ContourIndex(j, :-)

# 			coef = index(corr, i, j, b1=:-, b2=:-)
# 			t = __get_contour_op(lattice, ind1, ind2, op1, op2, coef)

# 			if !isnothing(mps2_mm)
# 				mps2_mm += apply!(t, vacuumstate(scalartype(hyb), lattice)) 
# 			else
# 				mps2_mm = apply!(t, vacuumstate(scalartype(hyb), lattice)) 
# 			end
# 		end
# 		# println(norm(mps_mm), " ", norm(mps2_mm), " ", dot(mps_mm, mps2_mm), " ", distance(mps_mm, mps2_mm))
# 		@test distance(mps_mm, mps2_mm) / norm(mps2_mm) < tol

# 		dt = 0.01

# 		mps_pp = dt * mps_pp + vacuumstate(lattice)
# 		mps_pm = dt * mps_pm + vacuumstate(lattice)
# 		mps_mp = dt * mps_mp + vacuumstate(lattice)
# 		mps_mm = dt * mps_mm + vacuumstate(lattice)


# 		mps2_pp, mps2_pm, mps2_mp, mps2_mm = influenceoperatorexponential(lattice, corr, dt, hyb, WII(), algexpan=algexpan)

# 		@test distance(mps_pp, mps2_pp) / norm(mps2_pp) < dt
# 		@test distance(mps_pm, mps2_pm) / norm(mps2_pm) < dt
# 		@test distance(mps_mp, mps2_mp) / norm(mps2_mp) < dt
# 		@test distance(mps_mm, mps2_mm) / norm(mps2_mm) < dt

# 		mps0 = mult!(mps2_pp, mps2_pm, trunc=trunc)
# 		mps0 = mult!(mps0, mps2_mp, trunc=trunc)
# 		mps0 = mult!(mps0, mps2_mm, trunc=trunc)

# 		for algmult in (SVDCompression(D=50), DMRGMult1(trunc=truncdimcutoff(D=50,ϵ=1.0e-6)))
# 			mps1 = differentialinfluencefunctional(lattice, corr, dt, hyb, WII(), algmult, algexpan=algexpan)
# 			_n = norm(mps1)
# 			mps2 = differentialinfluencefunctional(lattice, corr, dt, hyb, WI(), algmult, algexpan=algexpan)
# 			mps3 = differentialinfluencefunctional(lattice, corr, dt, hyb, ComplexStepper(WI()), algmult, algexpan=algexpan)
# 			mps4 = differentialinfluencefunctional(lattice, corr, dt, hyb, ComplexStepper(WII()), algmult, algexpan=algexpan)
# 			@test distance(mps1, mps0) / _n < dt
# 			@test distance(mps2, mps0) / _n < dt
# 			@test distance(mps3, mps0) / _n < dt
# 			@test distance(mps4, mps0) / _n < dt
# 		end

# 	end


# end
