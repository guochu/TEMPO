# MPO Hamiltonians, Schur MPO tensors, decay terms and steppers

function prodmpo(::Type{T}, ds::Vector{Int}, positions::Vector{Int}, ops::Vector{<:AbstractMatrix}) where {T <: Number}
	(length(positions) == length(ops)) || throw(DimensionMismatch("positions and ops size mismatch"))
	(length(Set(positions)) == length(positions)) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	L = length(ds)
	mpotensors = Vector{Array{T, 4}}(undef, L)
	for i in 1:L
		pos = findfirst(x->x==i, positions)
		dj = ds[i]
		if isnothing(pos)
			mj = one(zeros(T, dj, dj))
		else
			mj = ops[pos]
		end
		mpotensors[i] = reshape(mj, (1,dj,1,dj))
	end
	return ProcessTensor(mpotensors)
end
prodmpo(ds::Vector{Int}, positions::Vector{Int}, ops::Vector{<:AbstractMatrix}) = prodmpo(Float64, ds, positions, ops)
prodmpo(L::Int, positions::Vector{Int}, ops::Vector{<:AbstractMatrix}) = prodmpo([size(ops[1], 1) for _ in 1:L], positions, ops)

function longrange_xxz(J, Jzz, hz, α, p)
	sp, sm, z = p["+"], p["-"], p["z"]
	λ = exp(-α)
	return SchurMPOTensor(ExpDecayOpTerm(2*J*sp, sp', 1.0, λ), hz * z) +
		SchurMPOTensor(ExpDecayOpTerm(2*J*sm, sm', 1.0, λ)) +
		SchurMPOTensor(ExpDecayOpTerm(Jzz*z, z, 1.0, λ))
end

function longrange_xxz_ham(L, hz, J, Jzz, α, p)
	sp, sm, z = p["+"], p["-"], p["z"]
	mpo = prodmpo(L, [1], [hz * z])
	for i in 2:L
		mpo += prodmpo(L, [i], [hz * z])
	end
	canonicalize!(mpo)
	for i in 1:L
	    for j in i+1:L
	    	coeff = exp(-α*(j-i))
	    	mpo += prodmpo(L, [i, j], [2*J*coeff*sp, sp'])
	    	mpo += prodmpo(L, [i, j], [2*J*coeff*sm, sm'])
	    	mpo += prodmpo(L, [i, j], [Jzz*coeff*z, z])
	    	canonicalize!(mpo)
	    end
	end
	return mpo
end

longrange_xxz_mpoham(L, hz, J, Jzz, α, p) = MPOHamiltonian([longrange_xxz(J, Jzz, hz, α, p) for i in 1:L])

function powlaw_xxz(L, J, Jzz, hz, α, p)
	sp, sm, z = p["+"], p["-"], p["z"]
	return SchurMPOTensor(ExpDecayOpTerm(2*J*sp, sp', 1.0, exp(-1)), hz * z) +
		SchurMPOTensor(ExpDecayOpTerm(2*J*sm, sm', 1.0, exp(-1))) +
		SchurMPOTensor(expand_decayterm(PowerlawDecayTerm(z, Jzz*z, α=α), len=L, alg=OverDeterminedProny(tol=1.0e-8)))
end

powerlaw_xxz_mpoham(L, J, Jzz, hz, α, p) = MPOHamiltonian([powlaw_xxz(L, J, Jzz, hz, α, p) for i in 1:L])

function powerlaw_xxz_ham(L, J, Jzz, hz, α, p)
	sp, sm, z = p["+"], p["-"], p["z"]
	mpo = prodmpo(L, [1], [hz * z])
	for i in 2:L
		mpo += prodmpo(L, [i], [hz * z])
	end
	canonicalize!(mpo)
	for i in 1:L
	    for j in i+1:L
	    	coeff = exp(-(j-i))
	    	mpo += prodmpo(L, [i, j], [2*J*coeff*sp, sp'])
	    	mpo += prodmpo(L, [i, j], [2*J*coeff*sm, sm'])
	    	coeff = (j-i)^α
	    	mpo += prodmpo(L, [i, j], [Jzz*coeff*z, z])
	    	canonicalize!(mpo)
	    end
	end
	return mpo
end

@testset "decay terms                " begin
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
	# expansion packs into an ExpDecayOpSum (one channel per exponential)
	terms = expand_decayterm(g4)
	@test terms isa ExpDecayOpSum
	@test !isempty(terms.αs) && !isempty(terms.λs)
	@test length(terms.αs) == length(terms.λs)
	# function-valued decay requires a sampling length
	terms2 = expand_decayterm(g, len=10)
	@test !isempty(terms2.αs)
	# power-law convenience wrapper
	g5 = PowerlawDecayTerm(a, b; α=-2.0)
	@test g5 isa GenericDecayTerm
	# adjoint (vector-valued decay)
	@test scalartype(g4') == scalartype(g4)
end

@testset "MPOHamiltonian: long-range XXZ        " begin
	p = spin_half_matrices()
	for L in (2, 3)
		hz = 0.8
		J = 1
		Jzz = 1.2
		α = 0.9
		h1 = ProcessTensor(tompotensors(longrange_xxz_mpoham(L, hz, J, Jzz, α, p)))
		@test length(h1) == L
		h2 = longrange_xxz_ham(L, hz, J, Jzz, α, p)
		@test length(h2) == L
		@test distance(h1, h2) / norm(h2) < 1.0e-6
	end
end

@testset "MPOHamiltonian: power-law XXZ      " begin
	p = spin_half_matrices()
	L = 10
	α = -2.5
	hz = 0.8
	J = 1
	Jzz = 1.2
	h1 = ProcessTensor(tompotensors(powerlaw_xxz_mpoham(L, hz, J, Jzz, α, p)))
	h2 = powerlaw_xxz_ham(L, hz, J, Jzz, α, p)

	@test distance(h1, h2) / norm(h2) < 1.0e-5
end

@testset "MPOHamiltonian time evolution" begin
	p = spin_half_matrices()
	J, Jzz, hz, α = 1.0, 1.2, 0.8, 0.9
	m = longrange_xxz(J, Jzz, hz, α, p)
	h = MPOHamiltonian([m, m])
	dt = 1.0e-3
	h1 = timeevompo(h, dt, WI())
	h2 = timeevompo(h, dt, WII())
	h3 = timeevompo(h, dt)          # keyword version defaults to WII
	@test h1 isa MPO
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
	@test m1 isa Array            # single SchurMPOTensor evolves to a dense rank-4 tensor
	@test phydim(m1) == 2

	u1, u2 = timeevompo(h, dt, ComplexStepper(WI()))
	@test u1 isa MPO && u2 isa MPO
	d1, d2 = complex_stepper(dt)
	@test d1 + d2 ≈ dt
end
