# Cached observables: environments / expectationvalue / Zvalue / expectation /
# TransferMatrix / l_LL / r_RR on random states, plus rdm / quantummap / meanforcestate.
# Parametrized over ADT / PT and the three contours (two-point scans are sampled).

@testset "Cached observables: imaginary time" begin
	tol = 1.0e-6

	# ADT
	lattice = ADTLattice(N=3, δτ=0.1, contour=:imag)
	mps1 = randomadt(ComplexF64, length(lattice), D=4, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randomadt(ComplexF64, length(lattice), D=3, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mult(mps1, mps2)
	Zval = integrate(mps1, mps2)
	@test abs(integrate(mps) - Zval) / abs(Zval) < tol

	cache1 = environments(lattice, mps1, mps2)
	@test abs(Zvalue(cache1) - Zval) / abs(Zval) < tol
	cache2 = environments(lattice, mps)
	@test abs(Zvalue(cache2) - Zval) / abs(Zval) < tol

	op = randn(Float64, phydim(lattice))
	for i in 1:lattice.N
		t = ADTTerm(i, op)
		v1 = integrate(apply!(t, copy(mps))) / Zval
		@test abs(expectationvalue(t, cache1) - v1) / abs(v1) < tol
		@test abs(expectationvalue(t, cache2) - v1) / abs(v1) < tol
	end

	# PT
	lattice = PTLattice(N=3, δτ=0.1, contour=:imag)
	mps1 = randompt(ComplexF64, length(lattice), D=2, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randompt(ComplexF64, length(lattice), D=3, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mps1 * mps2
	Zval = integrate(lattice, mps1, mps2)
	@test abs(integrate(lattice, mps) - Zval) / abs(Zval) < tol

	cache1 = environments(lattice, mps1, mps2)
	@test abs(Zvalue(cache1) - Zval) / abs(Zval) < tol
	@test abs(Zvalue2(cache1) - Zval) / abs(Zval) < tol
	cache2 = environments(lattice, mps)
	@test abs(Zvalue(cache2) - Zval) / abs(Zval) < tol
	@test abs(Zvalue2(cache2) - Zval) / abs(Zval) < tol

	op = randn(Float64, phydim(lattice), phydim(lattice))
	for i in (1, 2)
		t = ProdFockTerm(i, op)
		v1 = integrate(lattice, apply!(t, copy(mps))) / Zval
		@test abs(expectationvalue(t, cache1) - v1) / abs(v1) < tol
		@test abs(expectationvalue(t, cache2) - v1) / abs(v1) < tol
	end
end

@testset "Cached observables: real time" begin
	tol = 1.0e-6

	# ADT
	lattice = ADTLattice(N=4, δt=0.1, contour=:real)
	mps1 = randomadt(ComplexF64, length(lattice), D=4, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randomadt(ComplexF64, length(lattice), D=3, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mult(mps1, mps2)
	Zval = integrate(mps1, mps2)
	cache1 = environments(lattice, mps1, mps2)
	cache2 = environments(lattice, mps)
	@test abs(Zvalue(cache1) - Zval) / abs(Zval) < tol
	@test abs(Zvalue(cache2) - Zval) / abs(Zval) < tol

	op1 = randn(Float64, phydim(lattice))
	op2 = randn(Float64, phydim(lattice))
	# sampled subset of two-point terms
	for i in (1, 2), b1 in branches(lattice), j in (1, 4), b2 in branches(lattice)
		pos1 = index(lattice, i, branch=b1)
		pos2 = index(lattice, j, branch=b2)
		(pos1 != pos2) || continue
		t = ADTTerm((pos1, pos2), (op1, op2))
		v1 = integrate(apply!(t, copy(mps))) / Zval
		@test abs(expectationvalue(t, cache1) - v1) / abs(v1) < tol
		@test abs(expectationvalue(t, cache2) - v1) / abs(v1) < tol
	end

	# PT (with a nontrivial initial state)
	lattice = PTLattice(N=4, δt=0.1, d=3, contour=:real)
	mps1 = randompt(ComplexF64, length(lattice), D=2, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randompt(ComplexF64, length(lattice), D=1, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mps1 * mps2
	ρ₀ = _rand_dm(phydim(lattice))
	mps2′ = initialstate!(copy(mps2), lattice, ρ₀)
	Zval = integrate(lattice, mps1, mps2′)
	mps′ = initialstate!(copy(mps), lattice, ρ₀)
	@test abs(integrate(lattice, mps′) - Zval) / abs(Zval) < tol

	cache1 = environments(lattice, mps1, mps2, ρ₀=ρ₀)
	@test abs(Zvalue(cache1) - Zval) / abs(Zval) < tol
	cache2 = environments(lattice, mps, ρ₀=ρ₀)
	@test abs(Zvalue(cache2) - Zval) / abs(Zval) < tol

	op1 = randn(Float64, phydim(lattice), phydim(lattice))
	op2 = randn(Float64, phydim(lattice), phydim(lattice))
	for i in (1, 2), b1 in branches(lattice), j in (1, 4), b2 in branches(lattice)
		pos1 = index(lattice, i, branch=b1)
		pos2 = index(lattice, j, branch=b2)
		(pos1 != pos2) || continue
		t = ProdFockTerm([pos1, pos2], [op1, op2])
		mps′ = apply!(t, copy(mps))
		mps′ = initialstate!(mps′, lattice, ρ₀)
		v1 = integrate(lattice, mps′) / Zval
		@test abs(expectationvalue(t, cache1) - v1) / abs(v1) < tol
		@test abs(expectationvalue(t, cache2) - v1) / abs(v1) < tol
	end
end

@testset "Cached observables: mixed time" begin
	tol = 1.0e-6

	# ADT
	lattice = ADTLattice(Nt=2, Nτ=3, δt=0.1, δτ=0.2, contour=:mixed)
	mps1 = randomadt(ComplexF64, length(lattice), D=4, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randomadt(ComplexF64, length(lattice), D=3, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mult(mps1, mps2)
	Zval = integrate(mps1, mps2)
	cache1 = environments(lattice, mps1, mps2)
	cache2 = environments(lattice, mps)
	@test abs(Zvalue(cache1) - Zval) / abs(Zval) < tol
	@test abs(Zvalue(cache2) - Zval) / abs(Zval) < tol

	op1 = randn(Float64, phydim(lattice))
	op2 = randn(Float64, phydim(lattice))
	for b1 in branches(lattice), i in (1, 2), b2 in branches(lattice)
		N2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
		pos1 = index(lattice, i, branch=b1)
		for j in (1, N2)
			pos2 = index(lattice, j, branch=b2)
			(pos1 != pos2) || continue
			t = ADTTerm((pos1, pos2), (op1, op2))
			v1 = integrate(apply!(t, copy(mps))) / Zval
			@test abs(expectationvalue(t, cache1) - v1) / abs(v1) < tol
			@test abs(expectationvalue(t, cache2) - v1) / abs(v1) < tol
		end
	end

	# PT
	lattice = PTLattice(Nt=3, Nτ=3, δt=0.1, δτ=0.2, d=2, contour=:mixed)
	mps1 = randompt(ComplexF64, length(lattice), D=2, d=phydim(lattice))
	canonicalize!(mps1)
	mps2 = randompt(ComplexF64, length(lattice), D=2, d=phydim(lattice))
	canonicalize!(mps2)
	mps = mps1 * mps2
	Zval = integrate(lattice, mps1, mps2)
	cache1 = environments(lattice, mps)
	cache2 = environments(lattice, mps1, mps2)
	@test abs(Zvalue(cache1) - Zval) / abs(Zval) < tol
	@test abs(Zvalue(cache2) - Zval) / abs(Zval) < tol

	op1 = randn(ComplexF64, phydim(lattice), phydim(lattice))
	op2 = randn(ComplexF64, phydim(lattice), phydim(lattice))
	for b1 in branches(lattice), i in (1, 2), b2 in branches(lattice)
		N2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
		pos1 = index(lattice, i, branch=b1)
		for j in (1, N2)
			pos2 = index(lattice, j, branch=b2)
			(pos1 != pos2) || continue
			t = ProdFockTerm([pos1, pos2], [op1, op2])
			v1 = integrate(lattice, apply!(t, copy(mps))) / Zval
			@test abs(expectationvalue(t, cache1) - v1) / abs(v1) < tol
			@test abs(expectationvalue(t, cache2) - v1) / abs(v1) < tol
		end
	end
end

@testset "expectation and TransferMatrix" begin
	N = 3
	tol = 1.0e-4
	lattice = ADTLattice(N=N, δt=0.1, contour=:real)
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
	@test m * r_RR(m) ≈ cache.hright[1]
	m12 = TransferMatrix(mps1, mps2)
	@test length(m12) == length(lattice)
	@test scaling(m12) ≈ scaling(mps1) * scaling(mps2)
	@test m12 * l_LL(m12) ≈ cache.hleft[end]
	m12j = TransferMatrix(2, mps1, mps2)
	@test length(m12j) == 1

	op = randn(Float64, 2)
	pos = index(lattice, 1, branch=:+)
	t = ADTTerm(pos, op)
	e1 = expectation(t, cache)
	@test e1 ≈ expectationvalue(t, cache) * Zvalue(cache)
	@test abs(e1 - integrate(apply!(t, copy(mps)))) / abs(e1) < tol
end

@testset "rdm, quantummap, meanforcestate" begin
	# rdm / quantummap on a real-time PT lattice
	lattice = PTLattice(N=3, δt=0.1, d=2, contour=:real)
	pt = randompt(ComplexF64, length(lattice), D=3, d=2)
	canonicalize!(pt)
	ρ = rdm(lattice, pt)
	@test size(ρ) == (2, 2)
	@test tr(ρ) ≈ integrate(lattice, pt)
	map1 = quantummap(lattice, pt)
	# the quantum map keeps the input/output indices untraced
	@test ndims(map1) == 4
	@test prod(size(map1)) == 16

	# mean-force state on an imaginary-time PT lattice
	lattice = PTLattice(N=3, δτ=0.1, contour=:imag)
	mps = randompt(ComplexF64, length(lattice), D=3, d=2)
	canonicalize!(mps)
	ρ = meanforcestate(lattice, mps)
	@test size(ρ) == (2, 2)
	@test ρ ≈ mfs(lattice, mps)
	@test tr(ρ) ≈ integrate(lattice, mps)
	mps2 = randompt(ComplexF64, length(lattice), D=2, d=2)
	canonicalize!(mps2)
	@test tr(meanforcestate(lattice, mps, mps2)) ≈ integrate(lattice, mps, mps2)
end
