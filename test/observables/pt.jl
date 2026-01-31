println("------------------------------------")
println("|           Process Tensor         |")
println("------------------------------------")

@testset "Cached observables: imaginary time" begin
	N = 3
	δτ = 0.1
	tol = 1.0e-6
	lattice = PTLattice(N=N, δτ=δτ, contour=:imag)
	mps1 = randompt(ComplexF64, length(lattice), D=2, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randompt(ComplexF64, length(lattice), D=3, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mps1 * mps2
	Zval = integrate(lattice, mps1, mps2)
	Zval1 = integrate(lattice, mps)
	@test abs(Zval1 - Zval) / abs(Zval) < tol

	
	cache1 = environments(lattice, mps1, mps2)
	Zval1 = Zvalue(cache1)
	@test abs(Zval1 - Zval) / abs(Zval) < tol
	Zval1 = Zvalue2(cache1)
	@test abs(Zval1 - Zval) / abs(Zval) < tol
	cache2 = environments(lattice, mps)
	Zval1 = Zvalue(cache2)
	@test abs(Zval1 - Zval) / abs(Zval) < tol
	Zval1 = Zvalue2(cache2)
	@test abs(Zval1 - Zval) / abs(Zval) < tol

	op = randn(Float64, phydim(lattice), phydim(lattice))
	for i in 1:1
		t = ProdFockTerm(i, op)
		mps′ = apply!(t, copy(mps))
		v1 = integrate(lattice, mps′) / Zval

		v2 = expectation(t, cache1)
		# println("v1=", v1, " v2=", v2)
		@test abs(v2 - v1) / abs(v1) < tol

		v2 = expectation(t, cache2)
		# println("v1=", v1, " v2=", v2)
		@test abs(v2 - v1) / abs(v1) < tol
	end

end

@testset "Cached diagonal observables: real time" begin
	N = 4
	δt = 0.1
	tol = 1.0e-6
	d = 3
	lattice = PTLattice(N=N, δt=δt, d=d, contour=:real)
	mps1 = randompt(ComplexF64, length(lattice), D=2, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randompt(ComplexF64, length(lattice), D=1, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mps1 * mps2
	ρ₀ = _rand_dm(phydim(lattice))
	mps2′ = initialstate!(copy(mps2), lattice, ρ₀)
	Zval = integrate(lattice, mps1, mps2′)
	mps′ = initialstate!(copy(mps), lattice, ρ₀)
	Zval1 = integrate(lattice, mps′)
	@test abs(Zval1 - Zval) / abs(Zval) < tol
	# canonicalize!(mps)

	
	cache1 = environments(lattice, mps1, mps2, ρ₀=ρ₀)
	Zval1 = Zvalue(cache1)
	# println("Zval=", Zval, " Zval1=", Zval1)
	@test abs(Zval1 * d - Zval) / abs(Zval) < tol

	cache2 = environments(lattice, mps, ρ₀=ρ₀)
	Zval1 = Zvalue(cache2)
	# println("Zval=", Zval, " Zval1=", Zval1)
	@test abs(Zval1 * d - Zval) / abs(Zval) < tol

	op1 = randn(Float64, phydim(lattice), phydim(lattice))
	op2 = randn(Float64, phydim(lattice), phydim(lattice))
	for i in 1:lattice.N, b1 in branches(lattice)
		pos1 = index(lattice, i, branch=b1)
		for j in 1:lattice.N, b2 in branches(lattice)
			pos2 = index(lattice, j, branch=b2)
			if pos1 != pos2
				t = ProdFockTerm([pos1, pos2], [op1, op2])
				mps′ = apply!(t, copy(mps))
				mps′ = initialstate!(mps′, lattice, ρ₀)
				v1 = integrate(lattice, mps′) / Zval
				v2 = expectation(t, cache1)
				@test abs(v2 - v1) / abs(v1) < tol
				v2 = expectation(t, cache2)
				@test abs(v2 - v1) / abs(v1) < tol
			end
		end
	end
end
