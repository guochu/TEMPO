println("------------------------------------")
println("|             PTLattice            |")
println("------------------------------------")


@testset "PTLattice: imaginary time M1N1" begin
	lattice = PTLattice(N=2, δτ=0.1, d=3, contour=:imag, ordering=M2M1())
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, ImagPTLattice)
	@test lattice.ordering == M2M1()
	@test scalartype(lattice) == Float64
	@test length(lattice) == 2
	@test lattice.d == 3
	@test lattice.N == 2
	@test lattice.β == 0.2
	@test lattice.δτ == 0.1
	@test lattice.τs == 0:0.1:0.2
	@test lattice.T == 5
	@test index(lattice, 2) == 1
	@test index(lattice, 1) == 2

	# mps = vacuumstate(lattice)
	# @test integrate(lattice, mps) ≈ 27 atol = 1.0e-6

end

@testset "PTLattice: real time M2m2M1m1" begin
	lattice = PTLattice(N=1, δt=0.1, contour=:real, ordering=M2m2M1m1())
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, RealPTLattice)
	@test lattice.ordering == M2m2M1m1()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 2
	@test lattice.d == 2
	@test lattice.N == 1
	@test lattice.t == 0.1
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 1, branch=:+) == 1
	@test index(lattice, 1, branch=:-) == 2

	# mps = vacuumstate(lattice)
	# @test integrate(lattice, mps) ≈ 16 atol = 1.0e-6

end

@testset "PTLattice: mixed time M2M1_m1M1m2M2" begin
	# one band
	lattice = PTLattice(Nt=1, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed, ordering=M2M1_m1M1m2M2())
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, MixedPTLattice)
	@test lattice.ordering == M2M1_m1M1m2M2()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 4
	@test lattice.Nt == 1
	@test lattice.Nτ == 2
	@test lattice.d == 2
	@test lattice.t == 0.05
	@test lattice.β == 0.2
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.05
	@test lattice.τs == 0:0.1:0.2

	# imaginary time axis
	@test index(lattice, 2, branch=:τ) == 1
	@test index(lattice, 1, branch=:τ) == 2

	# real time axis
	@test index(lattice, 1, branch=:-) == 3
	@test index(lattice, 1, branch=:+) == 4

	# mps = vacuumstate(lattice)
	# @test integrate(lattice, mps) ≈ 2^7 atol = 1.0e-6
end

