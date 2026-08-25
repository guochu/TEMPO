println("------------------------------------")
println("|               Basic              |")
println("------------------------------------")
println()

@testset "boundarycondition          " begin
	lattice = ADTLattice(N=3, δt=0.1, contour=:real)
	x = randomadt(ComplexF64, length(lattice), D=4, d=2)
	canonicalize!(x)
	xc = copy(x)
	ρimp = spin_up()
	x2 = boundarycondition(x, lattice, ρ₀=ρimp)
	@test x2 !== x
	@test distance(x2, boundarycondition!(copy(x), lattice, ρ₀=ρimp)) < 1.0e-12
	# the non-mutating version must not modify the input
	@test distance(x, xc) < 1.0e-14
end

@testset "sysdynamics!               " begin
	N = 4
	δt = 0.1
	lattice = ADTLattice(N=N, δt=δt, contour=:real)
	model = ImpurityHamiltonian(0.5 .* pauli_x())
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)
	a = sysdynamics!(vacuumstate(ComplexF64, lattice), lattice, model, trunc=trunc)
	# forward branch only differs from both-branch evolution
	b = sysdynamics!(vacuumstate(ComplexF64, lattice), lattice, model, trunc=trunc, branch=:+)
	@test distance(a, b) > 1.0e-6
	sysdynamics!(b, lattice, model, trunc=trunc, branch=:-)
	@test distance(a, b) < 1.0e-10
	# non-mutating wrapper consistency
	m1 = sysdynamics(lattice, model, trunc=trunc)
	m2 = sysdynamics!(vacuumstate(ComplexF64, lattice), lattice, model, trunc=trunc)
	@test distance(m1, m2) < 1.0e-12
end

@testset "hybriddynamics_naive!      " begin
	δt = 0.1
	N = 2
	β = 2
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β)
	lattice = ADTLattice(N=N, δt=δt, contour=:real)
	corr = correlationfunction(bath, lattice)
	hyb = AdditiveHyb(diag(pauli_z()))
	trunc = truncdimcutoff(D=30, ϵ=1.0e-12)
	m1 = hybriddynamics(lattice, corr, hyb, trunc=trunc)
	m2 = hybriddynamics_naive(lattice, corr, hyb, trunc=trunc)
	@test distance(m1, m2) / norm(m2) < 1.0e-3
	m3 = hybriddynamics_naive!(vacuumstate(ComplexF64, lattice), lattice, corr, hyb, trunc=trunc)
	@test distance(m2, m3) < 1.0e-12
end
