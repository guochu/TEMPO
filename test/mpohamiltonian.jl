println("------------------------------------")
println("|          MPO Hamiltonian         |")
println("------------------------------------")
println()

@testset "GenericDecayTerm           " begin
	a = randn(2, 2)
	b = randn(2, 2)
	f = x -> exp(-x)
	g = GenericDecayTerm(a, b, f)
	@test g isa GenericDecayTerm
	@test g.coeff == 1.0
	@test g.a == a && g.b == b
	g2 = GenericDecayTerm(a, b; f=f, coeff=2.0)
	@test g2.coeff == 2.0
	mid = one(2) .+ 0.1 .* randn(2, 2)
	g3 = GenericDecayTerm(a, b; middle=mid, f=f)
	@test g3.m == mid
	# vector-valued decay function
	fv = exp.(-(0:10))
	g4 = GenericDecayTerm(a, b, fv)
	@test g4 isa GenericDecayTerm
	terms = exponential_expansion(g4)
	@test !isempty(terms)
	@test all(t -> t isa ExponentialDecayTerm, terms)
	# function-valued decay requires a sampling length
	terms2 = exponential_expansion(g, len=10)
	@test !isempty(terms2)
	# power-law convenience wrapper
	g5 = PowerlawDecayTerm(a, b; α=-2.0)
	@test g5 isa GenericDecayTerm
	# adjoint (vector-valued decay)
	@test scalartype(g4') == scalartype(g4)
end

@testset "MPOHamiltonian time evolution" begin
	p = spin_half_matrices()
	sp, sm, z = p["+"], p["-"], p["z"]
	J, Jzz, hz, α = 1.0, 1.2, 0.8, 0.9
	m = SchurMPOTensor(hz*z, [ExponentialDecayTerm(2*J*sp, sp', α=exp(-α)),
	                          ExponentialDecayTerm(2*J*sm, sm', α=exp(-α)),
	                          ExponentialDecayTerm(Jzz*z, z, α=exp(-α))])
	h = MPOHamiltonian([m, m])
	dt = 1.0e-3
	h1 = timeevompo(h, dt, WI())
	h2 = timeevompo(h, dt, WII())
	h3 = timeevompo(h, dt)          # keyword version defaults to WII
	@test h1 isa MPOHamiltonian
	@test length(h1) == 2
	@test WI() isa FirstOrderStepper
	@test WII() isa FirstOrderStepper

	t1 = ProcessTensor(tompotensors(h1))
	t2 = ProcessTensor(tompotensors(h2))
	t3 = ProcessTensor(tompotensors(h3))
	# WI and WII agree to first order in dt
	@test distance(t1, t2) / norm(t1) < 1.0e-4
	@test distance(t2, t3) < 1.0e-12

	m1 = timeevompo(m, dt, WI())
	@test m1 isa SparseMPOTensor
	@test phydim(m1) == 2

	u1, u2 = timeevompo(h, dt, ComplexStepper(WI()))
	@test u1 isa MPOHamiltonian && u2 isa MPOHamiltonian
	d1, d2 = complex_stepper(dt)
	@test d1 + d2 ≈ dt
end
