# Dense/tensor linear algebra: orthogonal factorizations, tsvd, permute, isometry, distances

@testset "distance2                 " begin
	a = randn(3, 4)
	b = randn(3, 4)
	@test distance2(a, b) ≈ norm(a - b)^2
	@test distance(a, b) ≈ norm(a - b)
	@test distance2(a, a) == 0.0
	@test distance(a, a) == 0.0
end

@testset "isometry                  " begin
	@test isometry(3) == Matrix(I, 3, 3)
	@test isometry(3) isa Matrix{Float64}
	@test isometry(ComplexF64, 3) == Matrix(I, 3, 3)
	@test isometry(ComplexF64, 3) isa Matrix{ComplexF64}
	i34 = isometry(3, 4)
	@test i34 == [1 0 0 0; 0 1 0 0; 0 0 1 0]
	@test i34 * i34' == Matrix(I, 3, 3)
	i32 = isometry(ComplexF64, 3, 2)
	@test i32 == [1 0; 0 1; 0 0]
	@test i32' * i32 == Matrix(I, 2, 2)
end

@testset "permute                   " begin
	a = randn(3, 4)
	@test permute(a, (2, 1)) ≈ permutedims(a, (2, 1))
	@test permute(a, (1,), (2,)) ≈ a
	t = randn(2, 3, 4)
	@test permute(t, (1, 2), (3,)) ≈ permute(t, (1, 2, 3))
end

@testset "tsvd and tsvd!            " begin
	a = randn(ComplexF64, 6, 5)
	ac = copy(a)
	# in-place tensor version (also covers permute + reconstruction)
	t = randn(2, 3, 4)
	tc = copy(t)
	u, s, v, err = tsvd!(t, (1, 2), (3,))
	md = length(s)
	@test size(u) == (2, 3, md)
	@test size(v) == (md, 4)
	r = zeros(size(tc))
	for k in 1:md
		r += u[:, :, k] .* s[k] .* reshape(v[k, :], 1, 1, 4)
	end
	@test r ≈ tc
	@test err == 0.0

	for alg in (SVD(), SDD())
		u, s, v, err = tsvd(a; alg=alg)
		@test a == ac                       # input not modified
		@test err == 0.0
		@test u * Diagonal(s) * v ≈ a
		@test all(s .>= 0)
		# alg keyword also available for tsvd!
		u2, s2, v2, err2 = tsvd!(copy(a); alg=alg)
		@test s2 ≈ s
	end
	# SVD/SDD drivers give the same singular values
	_, s_svd, _, _ = tsvd(a; alg=SVD())
	_, s_sdd, _, _ = tsvd(a; alg=SDD())
	@test s_svd ≈ s_sdd
end

@testset "leftorth and rightorth    " begin
	A = randn(ComplexF64, 6, 4)
	Ac = copy(A)
	@test OrthogonalFactorizationAlgorithm <: Any
	# default algs
	Q, R = leftorth(A)
	@test A == Ac
	@test Q * R ≈ A
	@test Q' * Q ≈ Matrix(I, 4, 4)
	L, Q2 = rightorth(A)
	@test A == Ac
	@test L * Q2 ≈ A
	@test Q2 * Q2' ≈ Matrix(I, 4, 4)
	# all algorithms
	for alg in (QR(), QRpos(), SVD(), SDD(), Polar())
		Q, R = leftorth(A; alg=alg)
		@test A == Ac
		@test Q * R ≈ A
		@test Q' * Q ≈ Matrix(I, 4, 4)
	end
	for alg in (LQ(), LQpos(), SVD(), SDD())
		L, Q = rightorth(A; alg=alg)
		@test A == Ac
		@test L * Q ≈ A
		@test Q * Q' ≈ Matrix(I, 4, 4)
	end
	# Polar right-orthogonalization requires a wide matrix
	Aw = randn(ComplexF64, 4, 6)
	L, Q = rightorth(Aw; alg=Polar())
	@test L * Q ≈ Aw
	@test Q * Q' ≈ Matrix(I, 4, 4)
	# SVD-based orthogonalization with atol truncates small singular values
	B = Diagonal([1.0, 1.0e-12, 1.0])
	Qb, Rb = leftorth(Matrix(B); alg=SVD(), atol=1.0e-9)
	@test size(Qb, 2) == 2
	@test Qb * Rb ≈ B
	# tensor versions
	T = randn(4, 3, 2)
	Tc = copy(T)
	u, v = leftorth(T, (1, 2), (3,))
	@test T == Tc
	s = size(v, 1)
	@test size(u) == (4, 3, s)
	@test size(v) == (s, 2)
	@test reshape(u, :, s)' * reshape(u, :, s) ≈ Matrix(I, s, s)
	@test reshape(reshape(u, :, s) * reshape(v, s, :), 4, 3, 2) ≈ Tc
	T = randn(4, 3, 2)
	Tc = copy(T)
	u, v = rightorth(T, (1,), (2, 3))
	@test T == Tc
	s = size(v, 1)
	@test size(u) == (4, s)
	@test size(v) == (s, 3, 2)
	@test reshape(v, s, :) * reshape(v, s, :)' ≈ Matrix(I, s, s)
	@test reshape(reshape(u, :, s) * reshape(v, s, :), 4, 3, 2) ≈ Tc
end
