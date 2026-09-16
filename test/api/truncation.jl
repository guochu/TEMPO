# Truncation schemes, SVD compression and related utilities

@testset "truncation schemes          " begin
	@test truncdim(3) isa TruncateDim
	@test truncdim(D=4) isa TruncateDim
	@test truncdim(D=4).D == 4
	@test trunccutoff(ϵ=1.0e-3) isa TruncateCutoff
	@test trunccutoff(ϵ=1.0e-3).ϵ == 1.0e-3
	@test truncdimcutoff(D=5, ϵ=1.0e-3) isa TruncationDimCutoff
	@test truncdimcutoff(5, 1.0e-3) isa TruncationDimCutoff
	@test NoTruncation() isa TruncationScheme
	@test SVDCompression(truncdimcutoff(D=10, ϵ=1.0e-8)) isa SVDCompression

	a = randn(6, 5)
	# tsvd! destroys its input, use the non-mutating tsvd since a is reused below
	u, s, v, err = tsvd(a, trunc=truncdim(3))
	@test length(s) == 3
	@test err ≈ norm(svdvals(a)[4:end])
	u2, s2, v2, err2 = tsvd(a, trunc=trunccutoff(ϵ=1.0e-10))
	@test length(s2) == 5
	@test err2 < 1.0e-8
	u3, s3, v3, err3 = tsvd(a, trunc=truncdimcutoff(D=2, ϵ=1.0e-10))
	@test length(s3) == 2
	@test err3 > 0
	u4, s4, v4, err4 = tsvd(a)
	@test length(s4) == 5
	@test err4 == 0.0
	@test u4 * Diagonal(s4) * v4 ≈ a
	# tsvd! on a copy gives the same result
	u5, s5, v5, err5 = tsvd!(copy(a), trunc=truncdim(3))
	@test s5 ≈ s
end

@testset "renyi_entropy             " begin
	v = [0.25, 0.75]
	@test renyi_entropy(v) ≈ -(0.25 * log(0.25) + 0.75 * log(0.75))
	@test renyi_entropy(v; α=2) ≈ -log(0.25^2 + 0.75^2)
	@test renyi_entropy([1.0]) == 0.0
	@test_throws ArgumentError renyi_entropy([0.5, 0.6])     # not normalized
	@test_throws ArgumentError renyi_entropy([-0.5, 1.5])    # negative entries
	# on normalized squared singular values
	u, s, v2, _ = tsvd!(randn(5, 5))
	p = s.^2 ./ sum(s.^2)
	@test renyi_entropy(p) > 0
end
