println("------------------------------------")
println("|          MeanForceState          |")
println("------------------------------------")
println()

@testset "meanforcestate and mfs     " begin
	lattice = PTLattice(N=3, δτ=0.1, contour=:imag)
	mps = randompt(ComplexF64, length(lattice), D=3, d=2)
	canonicalize!(mps)
	ρ = meanforcestate(lattice, mps)
	@test size(ρ) == (2, 2)
	@test ρ ≈ mfs(lattice, mps)
	@test tr(ρ) ≈ integrate(lattice, mps)
	mps2 = randompt(ComplexF64, length(lattice), D=2, d=2)
	canonicalize!(mps2)
	@test tr(meanforcestate(lattice, mps, mps2)) ≈ integrate(lattice, mps, mps2)
end
