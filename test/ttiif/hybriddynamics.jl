println("------------------------------------")
println("|      TTI-IF Hybriddynamics       |")
println("------------------------------------")

@testset "TTI-IF-Hybriddynamics: imaginary-time" begin

	N = 6
	δτ = 0.1
	β = N * δτ
	rtol = 1.0e-2

	trunc = truncdimcutoff(D=50, ϵ=1.0e-6, add_back=0)

	base_alg = PartialIF(trunc=trunc)
	algs = [TranslationInvariantIF(k=5, fast=true), TranslationInvariantIF(k=5, fast=false)]

	spec = Leggett(d=1, ωc=1)

	bath = bosonicbath(spec, β=β, μ=0)
	corr = Δτ(bath, N=N, δτ=δτ)

	
	for d in (2,)
		lattice = ADTLattice(N=N, δτ=δτ, d=d, contour=:imag)

		hyb = AdditiveHyb(randn(Float64, d))

		mpsI1 = hybriddynamics(lattice, corr, hyb, base_alg) 

		for alg in algs
			mpsI2 = hybriddynamics(lattice, corr, hyb, alg)
			@test distance(mpsI1, mpsI2) / norm(mpsI1) < rtol
		end
	end
end


# @testset "TTI-IF-Hybriddynamics: real-time" begin

# 	N = 3
# 	δt = 0.1
# 	β = 1

# 	rtol = 1.0e-2
# 	trunc = truncdimcutoff(D=100, ϵ=1.0e-6, add_back=0)

# 	base_alg = PartialIF(trunc=trunc)
# 	alg2 = TranslationInvariantIF(k=5, algevo=WII(), algmult=SVDCompression(trunc))
# 	alg3 = TranslationInvariantIF(k=5, algmult=DMRGMult1(trunc=trunc, initguess=:svd))
# 	alg4 = TranslationInvariantIF(k=5, algmult=DMRGMult1(trunc=trunc, initguess=:pre))
# 	alg5 = TranslationInvariantIF(k=5, algmult=DMRGMult1(trunc=trunc, initguess=:rand, maxiter=10))
# 	alg6 = TranslationInvariantIF(k=5, algmult=DMRGMult1(trunc=trunc), fast=false)

# 	algs = [alg2, alg3, alg4, alg5, alg6]

# 	spec = Leggett(d=1, ωc=1)

# 	bath = bosonicbath(spec, β=β, μ=0)
# 	corr = Δt(bath, N=N, t=N*δt)

# 	d = 2

# 	# println("ordering is ", ordering)
# 	lattice = PTLattice(N=N, δt=δt, contour=:real, d=d)

# 	hyb = NonAdditiveHyb(_rand_ham(ComplexF64, d))

# 	mpsI1 = hybriddynamics_naive(lattice, corr, hyb, base_alg) 

# 	for alg in algs
# 		mpsI2 = hybriddynamics(lattice, corr, hyb, alg)
# 		@test distance(mpsI1, mpsI2) / norm(mpsI1) < rtol
# 	end


# end