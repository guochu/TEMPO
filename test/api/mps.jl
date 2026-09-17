# ADT and ProcessTensor MPS-level algebra, parametrized over the two tensor types
# to avoid duplicating the same checks (the former adt.jl / pt.jl testsets)

const MPSConstructors = ("ADT" => randomadt, "ProcessTensor" => randompt)

@testset "arithmetic and canonicalize" begin
	for (name, randmps) in MPSConstructors
		L = 6
		D = 6
		tol = 1.0e-7
		@testset "$name" begin
			for T in (Float64, ComplexF64)
				psi = randmps(T, L, D=D)
				@test scalartype(psi) == T
				@test space_l(psi) == 1
				@test space_r(psi) == 1

				@test bond_dimension(psi) <= D
				@test bond_dimensions(psi) isa Vector{Int}
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
	end
end

@testset "multiplications            " begin
	L = 6
	chi = 20
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	alg1 = SVDCompression(trunc)
	alg2 = DMRG1(trunc, initguess=:svd)
	alg3 = DMRG1(trunc, initguess=:rand, maxiter=10)
	alg4 = DMRG1(trunc, initguess=:pre, maxiter=10)
	# trunc is parameterized: `SVDCompression` accepts any TruncationScheme,
	# while `DMRG1` requires a scheme carrying the bond dimension `D`
	alg5 = SVDCompression(truncdim(chi))
	alg6 = DMRG1(truncdim(chi), initguess=:svd)
	alg7 = SVDCompression(truncrelerr(ϵ=1.0e-10))
	algs = [alg1, alg2, alg3, alg4, alg5, alg6, alg7]
	tol = 1.0e-7
	for (name, randmps) in MPSConstructors
		@testset "$name" begin
			for T in (Float64, ComplexF64)
				psi1 = randmps(T, L, D=4)
				psi2 = randmps(T, L, D=4)

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
	end
end

@testset "zip-up integration         " begin
	tol = 1.0e-7
	# ADT: zip-up without a lattice
	for T in (ComplexF64,)
		psi1 = randomadt(T, 6, D=4)
		psi2 = randomadt(T, 6, D=4)
		v1 = integrate(psi1 * psi2)
		v2 = integrate(psi1, psi2)
		@test abs(v2 - v1) / abs(v1) < tol
		canonicalize!(psi1)
		canonicalize!(psi2)
		v1 = integrate(psi1 * psi2)
		v2 = integrate(psi1, psi2)
		@test abs(v2 - v1) / abs(v1) < tol
	end
	# PT: zip-up on the three contours
	lattices = (PTLattice(N=4, δτ=0.1, d=2, contour=:imag),
	            PTLattice(N=4, δt=0.1, d=2, contour=:real),
	            PTLattice(Nt=2, δt=0.1, Nτ=4, δτ=0.1, d=2, contour=:mixed))
	for lattice in lattices
		L = length(lattice)
		pt1 = randompt(ComplexF64, L, d=2, D=3)
		canonicalize!(pt1)
		pt2 = randompt(ComplexF64, L, d=2, D=4)
		canonicalize!(pt2)
		v1 = integrate(lattice, pt1 * pt2)
		v2 = integrate(lattice, pt1, pt2)
		@test abs(v1 - v2) / abs(v1) < 1.0e-6
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

@testset "permute                     " begin
	# two sweeps of QR/SVD re-gauging introduce ~1e-8 numerical noise
	rtol = 1.0e-6
	# generous bond cutoff: the permutation itself is lossless
	trunc = truncdimcutoff(D=512, ϵ=1.0e-14, add_back=0)
	for (name, randmps) in MPSConstructors
		@testset "$name" begin
			for T in (Float64, ComplexF64), L in (4, 8)
				# mixed-canonical initial state
				psi = randmps(T, L, D=6)
				canonicalize!(psi, alg=Orthogonalize(trunc=trunc, normalize=false))
				@test iscanonical(psi)

				perm = randperm(L)
				psi1 = TEMPO.permute!(copy(psi), perm; trunc=trunc)
				# a) a permuted mixed-canonical state is still mixed-canonical
				@test iscanonical(psi1)
				# b) permuting again with the inverse permutation restores the original state
				@test distance(TEMPO.permute!(psi1, invperm(perm); trunc=trunc), psi) / norm(psi) < rtol

				# the identity permutation returns the state unchanged
				psi1 = TEMPO.permute(psi, collect(1:L); trunc=trunc)
				@test distance(psi1, psi) / norm(psi) < rtol
			end
		end
	end
end
