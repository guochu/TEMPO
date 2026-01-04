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
		lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)
		for T in (Float64, ComplexF64)
			hyb = NonAdditiveHyb(_rand_ham(T, d))

			mpsI1 = hybriddynamics_naive(lattice, corr, hyb, base_alg) 

			for alg in algs
				mpsI2 = hybriddynamics(lattice, corr, hyb, alg)
				@test distance(mpsI1, mpsI2) / norm(mpsI1) < rtol
			end

		end
	end

end