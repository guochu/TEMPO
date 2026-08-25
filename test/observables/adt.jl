println("------------------------------------")
println("|                 ADT              |")
println("------------------------------------")

@testset "Cached diagonal observables: imaginary time" begin
	N = 3
	δτ = 0.1
	tol = 1.0e-6
	lattice = ADTLattice(N=N, δτ=δτ, contour=:imag)
	mps1 = randomadt(ComplexF64, length(lattice), D=4, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randomadt(ComplexF64, length(lattice), D=3, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mult(mps1, mps2)
	Zval = integrate(mps1, mps2)
	Zval1 = integrate(mps)
	@test abs(Zval1 - Zval) / abs(Zval) < tol

	
	cache1 = environments(lattice, mps1, mps2)
	Zval1 = Zvalue(cache1)
	@test abs(Zval1 - Zval) / abs(Zval) < tol
	cache2 = environments(lattice, mps)
	Zval1 = Zvalue(cache2)
	@test abs(Zval1 - Zval) / abs(Zval) < tol

	op = randn(Float64, phydim(lattice))
	for i in 1:lattice.N
		t = ADTTerm(i, op)
		mps′ = apply!(t, copy(mps))
		v1 = integrate(mps′) / Zval
		v2 = expectationvalue(t, cache1)
		@test abs(v2 - v1) / abs(v1) < tol
		v2 = expectationvalue(t, cache2)
		@test abs(v2 - v1) / abs(v1) < tol
	end

end

@testset "Cached diagonal observables: real time" begin
	N = 4
	δt = 0.1
	tol = 1.0e-6
	lattice = ADTLattice(N=N, δt=δt, contour=:real)
	mps1 = randomadt(ComplexF64, length(lattice), D=4, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randomadt(ComplexF64, length(lattice), D=3, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mult(mps1, mps2)
	Zval = integrate(mps1, mps2)
	Zval1 = integrate(mps)
	@test abs(Zval1 - Zval) / abs(Zval) < tol

	
	cache1 = environments(lattice, mps1, mps2)
	Zval1 = Zvalue(cache1)
	@test abs(Zval1 - Zval) / abs(Zval) < tol
	cache2 = environments(lattice, mps)
	Zval1 = Zvalue(cache2)
	@test abs(Zval1 - Zval) / abs(Zval) < tol

	op1 = randn(Float64, phydim(lattice))
	op2 = randn(Float64, phydim(lattice))
	for i in 1:lattice.N, b1 in branches(lattice)
		pos1 = index(lattice, i, branch=b1)
		for j in 1:lattice.N, b2 in branches(lattice)
			pos2 = index(lattice, j, branch=b2)
			if pos1 != pos2
				t = ADTTerm((pos1, pos2), (op1, op2))
				mps′ = apply!(t, copy(mps))
				v1 = integrate(mps′) / Zval
				v2 = expectationvalue(t, cache1)
				@test abs(v2 - v1) / abs(v1) < tol
				v2 = expectationvalue(t, cache2)
				@test abs(v2 - v1) / abs(v1) < tol
			end
		end
	end
end

@testset "Cached diagonal observables: mixed time" begin
	Nτ = 3
	Nt = 2
	δt = 0.1
	tol = 1.0e-6
	lattice = ADTLattice(Nt=Nt, Nτ=Nτ, δt=δt, δτ=0.2, contour=:mixed)
	mps1 = randomadt(ComplexF64, length(lattice), D=4, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randomadt(ComplexF64, length(lattice), D=3, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mult(mps1, mps2)
	Zval = integrate(mps1, mps2)
	Zval1 = integrate(mps)
	@test abs(Zval1 - Zval) / abs(Zval) < tol

	
	cache1 = environments(lattice, mps1, mps2)
	Zval1 = Zvalue(cache1)
	@test abs(Zval1 - Zval) / abs(Zval) < tol
	cache2 = environments(lattice, mps)
	Zval1 = Zvalue(cache2)
	@test abs(Zval1 - Zval) / abs(Zval) < tol

	op1 = randn(Float64, phydim(lattice))
	op2 = randn(Float64, phydim(lattice))
	for b1 in branches(lattice)
		N1 = (b1==:τ) ? lattice.Nτ : lattice.Nt
		for i in 1:N1
			pos1 = index(lattice, i, branch=b1)
			for b2 in branches(lattice)
				N2 = (b2==:τ) ? lattice.Nτ : lattice.Nt
				for j in 1:N2
					pos2 = index(lattice, j, branch=b2)
					if pos1 != pos2
						t = ADTTerm((pos1, pos2), (op1, op2))
						mps′ = apply!(t, copy(mps))
						v1 = integrate(mps′) / Zval
						v2 = expectationvalue(t, cache1)
						@test abs(v2 - v1) / abs(v1) < tol
						v2 = expectationvalue(t, cache2)
						@test abs(v2 - v1) / abs(v1) < tol
					end					
				end
			end

		end

	end
end

@testset "expectation and TransferMatrix" begin
	N = 3
	δt = 0.1
	tol = 1.0e-4
	lattice = ADTLattice(N=N, δt=δt, contour=:real)
	mps1 = randomadt(ComplexF64, length(lattice), D=4, d=2)
	canonicalize!(mps1)
	mps2 = randomadt(ComplexF64, length(lattice), D=3, d=2)
	canonicalize!(mps2)
	mps = mult(mps1, mps2)
	Zval = integrate(mps)
	cache = environments(lattice, mps1, mps2)
	@test abs(Zvalue(cache) - Zval) / abs(Zval) < tol

	m = TransferMatrix(mps)
	@test length(m) == length(lattice)
	@test scaling(m) ≈ scaling(mps)
	@test only(l_LL(m) * m) ≈ integrate(mps)
	m12 = TransferMatrix(mps1, mps2)
	@test length(m12) == length(lattice)
	@test scaling(m12) ≈ scaling(mps1) * scaling(mps2)
	@test l_LL(m12) * m12 ≈ cache.hleft[end]
	m12j = TransferMatrix(2, mps1, mps2)
	@test length(m12j) == 1

	op = randn(Float64, 2)
	pos = index(lattice, 1, branch=:+)
	t = ADTTerm(pos, op)
	e1 = expectation(t, cache)
	@test e1 ≈ expectationvalue(t, cache) * Zvalue(cache)
	@test abs(e1 - integrate(apply!(t, copy(mps)))) / abs(e1) < tol
end