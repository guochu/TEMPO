println("------------------------------------")
println("|             FockLattice          |")
println("------------------------------------")


@testset "FockLattice: imaginary time M1N1" begin
	lattice = FockLattice(N=2, δτ=0.1, d=3, contour=:imag, ordering=M2M1())
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, ImagFockLattice)
	@test lattice.ordering == M2M1()
	@test scalartype(lattice) == Float64
	@test length(lattice) == 3
	@test lattice.d == 3
	@test lattice.N == 2
	@test lattice.k == 3
	@test lattice.β == 0.2
	@test lattice.δτ == 0.1
	@test lattice.τs == 0:0.1:0.2
	@test lattice.T == 5
	@test index(lattice, 3) == 1
	@test index(lattice, 2) == 2
	@test index(lattice, 1) == 3

	mps = vacuumstate(lattice)
	@test integrate(mps) ≈ 27 atol = 1.0e-6

end

@testset "FockLattice: real time M2m2M1m1" begin
	lattice = FockLattice(N=1, δt=0.1, contour=:real, ordering=M2m2M1m1())
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, RealFockLattice)
	@test lattice.ordering == M2m2M1m1()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 4
	@test lattice.d == 2
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 2, branch=:+) == 1
	@test index(lattice, 2, branch=:-) == 2
	@test index(lattice, 1, branch=:+) == 3
	@test index(lattice, 1, branch=:-) == 4

	mps = vacuumstate(lattice)
	@test integrate(mps) ≈ 16 atol = 1.0e-6

end

@testset "FockLattice: mixed time M2M1_m1M1m2M2" begin
	# one band
	lattice = FockLattice(Nt=1, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed, ordering=M2M1_m1M1m2M2())
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, MixedFockLattice)
	@test lattice.ordering == M2M1_m1M1m2M2()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 7
	@test lattice.Nt == 1
	@test lattice.kt == 2
	@test lattice.Nτ == 2
	@test lattice.kτ == 3
	@test lattice.d == 2
	@test lattice.t == 0.05
	@test lattice.β == 0.2
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.05
	@test lattice.τs == 0:0.1:0.2

	# imaginary time axis
	@test index(lattice, 3, branch=:τ) == 1
	@test index(lattice, 2, branch=:τ) == 2
	@test index(lattice, 1, branch=:τ) == 3

	# real time axis
	@test index(lattice, 1, branch=:-) == 4
	@test index(lattice, 1, branch=:+) == 5

	@test index(lattice, 2, branch=:-) == 6
	@test index(lattice, 2, branch=:+) == 7

	mps = vacuumstate(lattice)
	@test integrate(mps) ≈ 2^7 atol = 1.0e-6
end

