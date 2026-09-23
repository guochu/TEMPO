# SVDCompression / DMRG1 是 TEMPO 侧的算法配置（含 TEMPO 特有字段与构造）；
# 截断方案与张量分解（TruncationScheme、tsvd! 等）的测试由 FiniteMPSAlgorithms 负责

@testset "SVDCompression / DMRG1       " begin
	# `trunc` of `SVDCompression` accepts any TruncationScheme
	schemes = (truncdimcutoff(D=10, ϵ=1.0e-8), truncdim(10), truncrelerr(ϵ=1.0e-8), truncrelerr(1.0e-8), NoTruncation())
	for trunc in schemes
		@test SVDCompression(trunc) isa SVDCompression
		@test SVDCompression(trunc).trunc == trunc
	end
	# `trunc` of `DMRG1` must carry a maximum bond dimension `D`, which seeds the initial guess
	for trunc in (truncdimcutoff(D=10, ϵ=1.0e-8), truncdim(10))
		@test DMRG1(trunc) isa DMRG1
		@test DMRG1(trunc).trunc == trunc
	end
	@test_throws MethodError DMRG1(truncrelerr(ϵ=1.0e-8))
	@test_throws MethodError DMRG1(truncrelerr(1.0e-8))
	@test_throws MethodError DMRG1(NoTruncation())
	# keyword constructors and defaults
	@test SVDCompression() isa SVDCompression
	@test DMRG1() isa DMRG1
	@test SVDCompression(trunc=truncdim(4), verbosity=2).trunc == truncdim(4)
	@test DMRG1(trunc=truncdim(4), initguess=:rand).trunc == truncdim(4)
	# similar preserves the configuration and overrides single fields
	alg = SVDCompression(truncdim(4), verbosity=2)
	alg1 = similar(alg)
	@test alg1.trunc == truncdim(4) && alg1.verbosity == 2
	@test similar(alg; trunc=truncrelerr(ϵ=1.0e-8)).trunc isa TruncateRelError
	alg = DMRG1(truncdim(4), maxiter=7, initguess=:rand)
	alg1 = similar(alg)
	@test alg1.trunc == truncdim(4) && alg1.maxiter == 7 && alg1.initguess == :rand
	@test similar(alg; trunc=truncdimcutoff(D=6, ϵ=1.0e-8)).trunc isa TruncateDimCutoff
	# initguess validation
	@test_throws ArgumentError DMRG1(truncdim(4), initguess=:bad)
end
