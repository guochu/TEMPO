println("------------------------------------")
println("|           PT PartialIF           |")
println("------------------------------------")

@testset "PT PartialIF: imaginary time" begin
	δτ=0.1
	N = 2
	β = N * δτ
	tol = 1.0e-4

	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β)

	for d in [2,3]
		lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)
		corr = correlationfunction(bath, lattice)
		# complex Hermitian (real-symmetric case is covered in the real-time testset)
		op = (m = randn(ComplexF64, d, d); m + m')
		hyb = NonAdditiveHyb(op)
		for i in 1:lattice.Nτ
			ind1 = ContourIndex(i)
			p1 = partialif(lattice, ind1, corr, hyb)
			p2 = partialif_naive(lattice, ind1, corr, hyb)

			@test distance(p1, p2) / norm(p1) < tol
			# single-branch lattice: bond dimension d for any Hermitian operator
			@test maximum(bond_dimensions(p1)) == d
		end
	end

end


@testset "PT PartialIF: real time" begin
	δt=0.1
	N = 3
	t = N * δt
	β = 2
	tol = 1.0e-4

	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β)

	for d in [2,3]
		lattice = PTLattice(N=N, δt=δt, d=d, contour=:real)
		corr = correlationfunction(bath, lattice)
		# symmetric operator: bond dimension d
		op = (m = randn(d, d); m + m')
		hyb = NonAdditiveHyb(op)
		for i in 1:lattice.Nt, b1 in branches(lattice)
			ind1 = ContourIndex(i, branch=b1)
			p1 = partialif(lattice, ind1, corr, hyb)
			p2 = partialif_naive(lattice, ind1, corr, hyb)

			@test distance(p1, p2) / norm(p1) < tol
			@test maximum(bond_dimensions(p1)) == d
		end
		# non-symmetric Hermitian operator: still bond dimension d
		opc = (m = randn(ComplexF64, d, d); m + m')
		hybc = NonAdditiveHyb(opc)
		for i in 1:lattice.Nt, b1 in branches(lattice)
			ind1 = ContourIndex(i, branch=b1)
			p1 = partialif(lattice, ind1, corr, hybc)
			p2 = partialif_naive(lattice, ind1, corr, hybc)

			@test distance(p1, p2) / norm(p1) < tol
			@test maximum(bond_dimensions(p1)) == d
		end
	end
end


@testset "PT PartialIF: mixed time" begin
	δτ=0.1
	Nτ = 3
	δt=0.05
	Nt = 2
	t = Nt * δt
	β = Nτ * δτ
	tol = 1.0e-4

	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β)

	for d in [2,3]
		lattice = PTLattice(Nτ=Nτ, δτ=δτ, Nt=Nt, δt=δt, d=d, contour=:mixed)
		corr = correlationfunction(bath, lattice)
		# symmetric operator: bond dimension d
		op = (m = randn(d, d); m + m')
		hyb = NonAdditiveHyb(op)
		for b1 in branches(lattice)
			N = (b1 == :τ) ? lattice.Nτ : lattice.Nt
			for i in 1:N
				ind1 = ContourIndex(i, branch=b1)
				p1 = partialif(lattice, ind1, corr, hyb)
				p2 = partialif_naive(lattice, ind1, corr, hyb)

				@test distance(p1, p2) / norm(p1) < tol
				@test maximum(bond_dimensions(p1)) == d
			end
		end
		# non-symmetric Hermitian operator: still bond dimension d
		opc = (m = randn(ComplexF64, d, d); m + m')
		hybc = NonAdditiveHyb(opc)
		for b1 in branches(lattice)
			N = (b1 == :τ) ? lattice.Nτ : lattice.Nt
			for i in 1:N
				ind1 = ContourIndex(i, branch=b1)
				p1 = partialif(lattice, ind1, corr, hybc)
				p2 = partialif_naive(lattice, ind1, corr, hybc)

				@test distance(p1, p2) / norm(p1) < tol
				@test maximum(bond_dimensions(p1)) == d
			end
		end
	end
end


@testset "PT PartialIF: influence functional" begin
	δτ=0.1
	N = 2
	β = N * δτ
	tol = 1.0e-4

	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β)
	lattice = PTLattice(N=N, δτ=δτ, d=2, contour=:imag)
	corr = correlationfunction(bath, lattice)
	op = (m = randn(2, 2); m + m')
	hyb = NonAdditiveHyb(op)

	if1 = hybriddynamics(lattice, corr, hyb, PartialIF())
	if2 = hybriddynamics_naive(lattice, corr, hyb, PartialIF())

	@test distance(if1, if2) / norm(if1) < tol
end
