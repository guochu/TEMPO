println("------------------------------------")
println("|                ADT               |")
println("------------------------------------")

@testset "ADT: arithmetic and canonicalize" begin
	L = 6
	D = 6
	tol = 1.0e-7
	for T in (Float64, ComplexF64)
		psi = randomadt(T, L, D=D)
		@test scalartype(psi) == T
		@test space_l(psi) == 1
		@test space_r(psi) == 1

		@test bond_dimension(psi) <= D
		psi1 = leftorth!(deepcopy(psi), alg = Orthogonalize(QR(), normalize=false))
		@test norm(psi) ≈ norm(psi1) atol = tol
		@test distance(psi, psi1) / norm(psi) < tol

		psi1 = rightorth!(deepcopy(psi), alg = Orthogonalize(QR(), normalize=false))
		@test norm(psi) ≈ norm(psi1) atol = tol
		@test distance(psi, psi1) / norm(psi) < tol

		psi1 = leftorth!(deepcopy(psi), alg = Orthogonalize(QR(), normalize=true))
		@test isleftcanonical(psi1)
		psi1 = rightorth!(deepcopy(psi), alg = Orthogonalize(SVD(), normalize=true))
		@test isrightcanonical(psi1)
		psi1 = canonicalize!(deepcopy(psi), alg = Orthogonalize(SVD(), normalize=true))
		@test iscanonical(psi1)
		@test norm(2 * psi1) ≈ 2
		@test norm(psi1 / 2) ≈ 0.5
		@test norm(psi1 - psi1) ≈ 0. atol = tol
		@test distance(psi, psi) ≈ 0. atol = tol

		psi1 = canonicalize!(deepcopy(psi), alg=Orthogonalize(trunc=NoTruncation(), normalize=false))
		@test norm(psi) ≈ norm(psi1) atol = tol
		@test distance(psi, psi1) / norm(psi) < tol

	end
	
end


@testset "ADT: multiplications" begin
	L = 6
	chi = 20
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	alg1 = SVDCompression(trunc)
	alg2 = DMRGMult1(trunc, initguess=:svd)
	alg3 = DMRGMult1(trunc, initguess=:rand, maxiter=10)
	alg4 = DMRGMult1(trunc, initguess=:pre, maxiter=10)
	algs = [alg1, alg2, alg3, alg4]
	tol = 1.0e-7
	for T in (Float64, ComplexF64)
		psi1 = randomadt(T, L, D=4)
		psi2 = randomadt(T, L, D=4)

		psi3 = psi1 * psi2
		_n = norm(psi3)
		psi4 = mult(psi1, psi2, trunc=trunc)
		@test distance(psi3, psi4) / _n < tol


		canonicalize!(psi1)
		canonicalize!(psi2)
		psi5 = psi1 * psi2
		@test distance(psi3, psi5) / _n < tol

		psi4 = mult(psi1, psi2, trunc=trunc)
		@test distance(psi3, psi4) / _n < tol

		for alg in algs
			psi4 = mult(psi1, psi2, alg)
			@test distance(psi3, psi4) / _n < tol
		end
	end
end


@testset "ADT: zipup integrate" begin

	L = 6
	tol = 1.0e-7
	for T in (Float64, ComplexF64)
		psi1 = randomadt(T, L, D=4)
		psi2 = randomadt(T, L, D=4)

		psi3 = psi1 * psi2

		v1 = integrate(psi3)
		v2 = integrate(psi1, psi2)
		@test abs(v2 - v1) / abs(v1) < tol

		canonicalize!(psi1)
		canonicalize!(psi2)
		psi3 = psi1 * psi2

		v1 = integrate(psi3)
		v2 = integrate(psi1, psi2)
		@test abs(v2 - v1) / abs(v1) < tol
		
	end

end

@testset "accessors: scaling, phydims, indexmappings" begin
	lattice = ADTLattice(N=3, δt=0.1, contour=:real)
	mps = randomadt(ComplexF64, length(lattice), D=4, d=2)
	@test scaling(mps) isa Float64
	@test all(phydims(mps) .== 2)
	@test phydim(mps, 1) == 2
	@test phydim(mps, length(lattice)) == 2

	map_ = indexmappings(lattice)
	@test length(map_) == 2 * lattice.k
	for i in 1:lattice.k, f in (:+, :-)
		@test map_[(i, f)] == index(lattice, i, branch=f)
	end

	latt2 = ADTLattice(N=3, δτ=0.1, contour=:imag)
	map2 = indexmappings(latt2)
	@test length(map2) == latt2.k
	for i in 1:latt2.k
		@test map2[(i, :τ)] == index(latt2, i)
	end

	mps2 = randomadt(ComplexF64, length(lattice), D=3, d=2)
	canonicalize!(mps2)
	@test scaling(mps, mps2) ≈ scaling(mps) * scaling(mps2)
	m = mult(mps, mps2)
	# the total contraction is invariant under the renormalization performed inside mult!
	@test abs(integrate(m) - integrate(mps, mps2)) / abs(integrate(mps, mps2)) < 1.0e-4
	@test scaling(TransferMatrix(mps, mps2)) ≈ scaling(mps) * scaling(mps2)
end
