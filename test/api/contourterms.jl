# Contour operators and influence-functional terms: ADTTerm, FockTerm(s), apply!, hybridizations

@testset "ADTTerm / apply!" begin
	lattice = ADTLattice(N=3, δt=0.1, contour=:real)
	mps = randomadt(ComplexF64, length(lattice), D=4, d=2)
	canonicalize!(mps)
	Zval = integrate(mps)

	# single-site diagonal term (vector form) on a lattice position
	z = [0.6, -0.4]
	pos = index(lattice, 2, branch=:+)
	t = ADTTerm((pos,), z)
	@test t isa ADTTerm
	v1 = integrate(apply!(t, copy(mps))) / Zval
	# inserting the identity at a second position must not change the value
	pos2 = index(lattice, 1, branch=:+)
	v2 = integrate(apply!(ADTTerm((pos, pos2), (z, [1.0, 1.0])), copy(mps))) / Zval
	@test v1 ≈ v2 rtol = 1.0e-12
	# two-site term (tuple of diagonal vectors)
	z2 = [-0.3, 0.8]
        pos3 = index(lattice, 1, branch=:-)
        t2 = ADTTerm((pos, pos3), (z, z2))
        @test t2 isa ADTTerm
        v3 = integrate(apply!(t2, copy(mps))) / Zval
        # the constructor sorts the positions and reorders the operators accordingly
        t2b = ADTTerm((pos3, pos), (z2, z))
        v3b = integrate(apply!(t2b, copy(mps))) / Zval
        @test v3 ≈ v3b rtol = 1.0e-12
end

@testset "hybridizations and pairop" begin
	z = [-1.0 0.0; 0.0 1.0]
	h1 = AdditiveHyb(randn(2))
	h2 = NonAdditiveHyb(z)
	h3 = NonDiagonalHyb(randn(ComplexF64, 2, 2))
	@test h1 isa HybridizationStyle
	@test h2 isa HybridizationStyle
	@test h3 isa HybridizationStyle
	# AdditiveHyb diagonalizes into a pair of identical diagonal operators
	op1, op2 = pairop(h1)
	@test op1 == op2
	@test Diagonal(h1.op) ≈ op1
	op1, op2 = pairop(h2)
	@test op1 === h2.op && op2 === h2.op
	op1, op2 = pairop(h3)
	@test op2 == op1'
end

@testset "FockTerm / FockTermS / ProdFockTerm" begin
	lattice = PTLattice(N=3, δτ=0.1, contour=:imag)
	mps = randompt(ComplexF64, length(lattice), D=3, d=2)
	canonicalize!(mps)
	cache = environments(lattice, mps)
	Zval = Zvalue(cache)

	z = pauli_z()
	pos = 1
	ft = FockTermS(pos, z)
	@test ft isa FockTermS
	@test ft isa AbstractFockTerm
	ft3 = FockTerm([pos], [reshape(ComplexF64.(z), 1, 2, 1, 2)])
	@test ft3 isa FockTerm
	# FockTermS and FockTerm share the same application machinery
	v1 = expectationvalue(ft, cache)
	v2 = expectationvalue(ft3, cache)
	@test v1 ≈ v2
	# self-consistent reference: apply the term and integrate
	ref = integrate(lattice, apply!(ft, copy(mps))) / Zval
	@test abs(v1 - ref) / abs(ref) < 1.0e-6
	# a ContourOperator maps to a ProdFockTerm at the same lattice position
	c = ContourIndex(1, branch=:τ)
	p1 = lattice[c]
	@test expectationvalue(ContourOperator(c, z), cache) ≈ expectationvalue(ProdFockTerm(p1, z), cache)
end
