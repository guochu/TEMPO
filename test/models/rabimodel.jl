println("------------------------------------")
println("|            Rabi Model            |")
println("------------------------------------")

function rabi_ham(Ω; d)
	x = pauli_x()
	z = pauli_z()
	# sp = Array{Float64, 2}([0 0; 1 0])
	# ns = [0 0; 0 1.]
	Is = one(x)
	a = bosonaoperator(d=d)
	adag = a'
	n = bosondensityoperator(d=d)
	Ib = one(n)

	Himp = Ω * kron(x, Ib)
	Hbath = kron(Is, n)
	Hhyb = kron(z, adag+a)

	H = Himp + Hhyb + Hbath
	return H, n
end

@testset "Rabi model: imaginary-time" begin

	Ω = 0.5
	N = 20
	δτ = 0.1
	β = N * δτ
	chi = 100
	d = 50
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(N = N, δτ=δτ, contour=:imag)

	x = [0 1; 1 0]
	xop = Ω .* x
	z = [-1 0; 0 1]
	Is = one(x)
	Ib = one(zeros(d, d))
	model = BosonicImpurity(xop)

	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, trunc=trunc)
	
	bs = AdditiveHyb([z[i,i] for i in 1:size(z,1)])

	spec = DiracDelta(1)

	bath = bosonicbath(spec, β=β)

	corr = correlationfunction(bath, lattice)

	mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI′ = hybriddynamics_naive(lattice, corr, bs, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI′) < tol
	mps = mult!(mpsK, mpsI, trunc=trunc)


	H, Hbarebath = rabi_ham(Ω, d=d)

	ρ = exp(-β * H)
	z1 = integrate(mps)
	@test abs(z1 - tr(ρ) / tr(exp(-β .* Hbarebath))) / abs(z1) < tol


	## diagonal observables
	op = [-0.73 0; 0 0.5]
	zdiag = [op[i,i] for i in 1:size(z, 1)]

	pos1 = index(lattice, 1)
	t = ADTTerm(pos1, zdiag .* zdiag )
	mps2 = apply!(t, copy(mps))
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	for i in 2:N
		pos2 = index(lattice, i)
		t = ADTTerm((pos2,pos1), (zdiag, zdiag))
		# t = ADTTerm((i,1), reshape(kron(zdiag, zdiag), 2, 2))
		mps2 = apply!(t, copy(mps))
		v = integrate(mps2) / integrate(mps)
		push!(corrs, v)
	end
	

	A = kron(op, Ib)

	corrs2 = correlation_2op_1τ(H, A, A, 0:δτ:β, β=β)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	## off-diagonal observables
	op1 = [0 0; 0.7 0]
	op2 = [0 0.8;0 0 ]

	c1 = ContourIndex(1)

	ct = ContourOperator(c1, op1 * op2)
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, trunc=trunc)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	for i in 2:N
		c2 = ContourIndex(i)
		ct = ContourOperator([c2, c1], [op2, op1])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, trunc=trunc)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / integrate(mps)
		push!(corrs, v)
	end
	

	A1 = kron(op1, Ib)
	A2 = kron(op2, Ib)

	corrs2 = correlation_2op_1τ(H, A1, A2, 0:δτ:β, β=β)
	corrs2 = corrs2[1:length(corrs)]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end