# partialif vs partialif_naive on the three contours, parametrized over ADT / PT lattices

@testset "PartialIF vs naive: ADT lattice" begin
	spec = Leggett(d=1, ωc=1)
	tol = 1.0e-4

	# imaginary time
	δτ = 0.1
	N = 5
	bath = bosonicbath(spec, β=N*δτ)
	corr = Δτ(bath, N=N, δτ=δτ)
	for d in (2, 3)
		lattice = ADTLattice(N=N, δτ=δτ, d=d, contour=:imag)
		op = rand(d)
		hyb = AdditiveHyb(op)
		for i in 1:lattice.N
			ind1 = ContourIndex(i)
			p1 = partialif(lattice, ind1, corr, hyb)
			p2 = partialif_naive(lattice, ind1, corr, hyb)
			@test distance(p1, p2) / norm(p1) < tol
		end
	end

	# real time
	δt = 0.1
	N = 5
	bath = bosonicbath(spec, β=2)
	corr = Δt(bath, N=N, t=N*δt)
	d = 2
	lattice = ADTLattice(N=N, δt=δt, d=d, contour=:real)
	hyb = AdditiveHyb(rand(d))
	for i in 1:lattice.N, b1 in branches(lattice)
		ind1 = ContourIndex(i, branch=b1)
		p1 = partialif(lattice, ind1, corr, hyb)
		p2 = partialif_naive(lattice, ind1, corr, hyb)
		@test distance(p1, p2) / norm(p1) < tol
	end

	# mixed time
	δτ = 0.1
	Nτ = 5
	δt = 0.05
	Nt = 3
	bath = bosonicbath(spec, β=2)
	corr = Δm(bath, Nτ=Nτ, Nt=Nt, t=Nt*δt)
	d = 2
	lattice = ADTLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, d=d, contour=:mixed)
	hyb = AdditiveHyb(rand(d))
	for b1 in branches(lattice)
		M = (b1 == :τ) ? lattice.Nτ : lattice.Nt
		for i in 1:M
			ind1 = ContourIndex(i, branch=b1)
			p1 = partialif(lattice, ind1, corr, hyb)
			p2 = partialif_naive(lattice, ind1, corr, hyb)
			@test distance(p1, p2) / norm(p1) < tol
		end
	end
end

@testset "PartialIF vs naive: PT lattice" begin
	spec = Leggett(d=1, ωc=1)
	tol = 1.0e-4

	# imaginary time (complex Hermitian operator)
	δτ = 0.1
	N = 2
	bath = bosonicbath(spec, β=N*δτ)
	for d in (2, 3)
		lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)
		corr = correlationfunction(bath, lattice)
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

	# real time (symmetric and non-symmetric Hermitian operators)
	δt = 0.1
	N = 3
	bath = bosonicbath(spec, β=2)
	for d in (2, 3)
		lattice = PTLattice(N=N, δt=δt, d=d, contour=:real)
		corr = correlationfunction(bath, lattice)
		for op in ((m = randn(d, d); m + m'), (m = randn(ComplexF64, d, d); m + m'))
			hyb = NonAdditiveHyb(op)
			for i in 1:lattice.Nt, b1 in branches(lattice)
				ind1 = ContourIndex(i, branch=b1)
				p1 = partialif(lattice, ind1, corr, hyb)
				p2 = partialif_naive(lattice, ind1, corr, hyb)
				@test distance(p1, p2) / norm(p1) < tol
				@test maximum(bond_dimensions(p1)) == d
			end
		end
	end

	# mixed time (symmetric and non-symmetric Hermitian operators)
	δτ = 0.1
	Nτ = 3
	δt = 0.05
	Nt = 2
	bath = bosonicbath(spec, β=Nτ*δτ)
	d = 2
	lattice = PTLattice(Nτ=Nτ, δτ=δτ, Nt=Nt, δt=δt, d=d, contour=:mixed)
	corr = correlationfunction(bath, lattice)
	for op in ((m = randn(d, d); m + m'), (m = randn(ComplexF64, d, d); m + m'))
		hyb = NonAdditiveHyb(op)
		for b1 in branches(lattice)
			M = (b1 == :τ) ? lattice.Nτ : lattice.Nt
			for i in 1:M
				ind1 = ContourIndex(i, branch=b1)
				p1 = partialif(lattice, ind1, corr, hyb)
				p2 = partialif_naive(lattice, ind1, corr, hyb)
				@test distance(p1, p2) / norm(p1) < tol
				@test maximum(bond_dimensions(p1)) == d
			end
		end
	end
end
