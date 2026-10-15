# SVDCompression / DMRG1 是 TEMPO 侧的算法配置（含 TEMPO 特有字段与构造）；
# 截断方案与张量分解（TruncationScheme、tsvd! 等）的测试由 FiniteMPSAlgorithms 负责

@testset "SVDCompression / DMRG1       " begin
	# `trunc` of `SVDCompression` accepts any TruncationScheme
	schemes = (truncdimcutoff(D=10, ϵ=1.0e-8), truncdim(10), truncrelerr(ϵ=1.0e-8), truncrelerr(1.0e-8), NoTruncation())
	for trunc in schemes
		@test SVDCompression(trunc=trunc) isa SVDCompression
		@test SVDCompression(trunc=trunc).trunc == trunc
	end
	# `trunc` of `DMRG1` must carry a maximum bond dimension `D`, which seeds the initial guess
	for trunc in (truncdimcutoff(D=10, ϵ=1.0e-8), truncdim(10))
		@test DMRG1(trunc=trunc) isa DMRG1
		@test DMRG1(trunc=trunc).trunc == trunc
	end
	@test_throws MethodError DMRG1(trunc=truncrelerr(ϵ=1.0e-8))
	@test_throws MethodError DMRG1(trunc=truncrelerr(1.0e-8))
	@test_throws MethodError DMRG1(trunc=NoTruncation())
	# keyword constructors and defaults
	@test SVDCompression() isa SVDCompression
	@test DMRG1() isa DMRG1
	@test SVDCompression(trunc=truncdim(4), verbosity=2).trunc == truncdim(4)
	@test DMRG1(trunc=truncdim(4), initguess=:rand).trunc == truncdim(4)
	# `truncrelerr`/`NoTruncation` are rejected by both DMRG1 and DMRG2
	@test_throws MethodError DMRG2(trunc=truncrelerr(ϵ=1.0e-8))
	@test_throws MethodError DMRG2(trunc=NoTruncation())
	# DMRG2 keyword constructor and defaults (initguess: :svd or :pre)
	@test DMRG2() isa DMRG2
	@test DMRG2(trunc=truncdim(4), maxiter=9).trunc == truncdim(4)
	@test DMRG2(trunc=truncdimcutoff(D=10, ϵ=1.0e-8)).trunc == truncdimcutoff(D=10, ϵ=1.0e-8)
	@test DMRG2(trunc=truncdim(4), initguess=:pre).initguess == :pre
end
