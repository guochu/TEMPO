println("------------------------------------")
println("|            Rabi Model            |")
println("------------------------------------")


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
	model = ImpurityHamiltonian(xop)

	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	
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
	ids2N = [k for k in sampleidx(N) if k > 1]
	idsallN = [1; ids2N]
	for i in ids2N
		pos2 = index(lattice, i)
		t = ADTTerm((pos2,pos1), (zdiag, zdiag))
		# t = ADTTerm((i,1), reshape(kron(zdiag, zdiag), 2, 2))
		mps2 = apply!(t, copy(mps))
		v = integrate(mps2) / integrate(mps)
		push!(corrs, v)
	end
	

	A = kron(op, Ib)

	corrs2 = correlation_2op_1τ(H, A, A, 0:δτ:β, β=β)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	## off-diagonal observables
	op1 = [0 0; 0.7 0]
	op2 = [0 0.8;0 0 ]

	c1 = ContourIndex(1)

	ct = ContourOperator(c1, op1 * op2)
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	c2 = ContourIndex(1)
	for i in ids2N
		c1 = ContourIndex(i)
		ct = ContourOperator([c1, c2], [op1, op2])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / integrate(mps)
		push!(corrs, v)
	end
	

	A1 = kron(op1, Ib)
	A2 = kron(op2, Ib)

	corrs2 = correlation_2op_1τ(H, A1, A2, 0:δτ:β, β=β)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end

@testset "Rabi model: real-time" begin

	Ω = 0.5
	N = 10
	δt = 0.05
	β = 2
	t = N * δt
	chi = 100
	d = 50
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(N = N, δt=δt, contour=:real)

	# x = [0 1; 1 0]
	x = Matrix{ComplexF64}([0 im; -im 0])
	hop = Ω .* x
	z = [-1 0; 0 1]
	Is = one(x)
	Ib = one(zeros(d, d))
	model = ImpurityHamiltonian(hop)

	Hbarebath = bosondensityoperator(d=d)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib) + kron(Is, Hbarebath) + kron(z, a' + a)


	bs = AdditiveHyb([z[i,i] for i in 1:size(z,1)])
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI′ = hybriddynamics_naive(lattice, corr, bs, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI′) < tol


	ρ1 = zeros(2,2)
	ρ1[1,1] = 1
	ρ2 = 0.5 .* one(ρ1)

	for ρimp in [ρ1,ρ2]

		mpsK = sysdynamics(lattice, model, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
		mps = mult!(mpsK, mpsI, trunc=trunc)

		
		ρ = kron(ρimp, exp(-β * Hbarebath)) 

		## diagonal observables
		op = [-0.73 0; 0 0.5]
		zdiag = [op[i,i] for i in 1:size(z, 1)]

		pos1 = index(lattice, 1, branch=:+)
		m = ADTTerm(pos1, zdiag .* zdiag )
		mps2 = apply!(m, copy(mps))
		v = integrate(mps2) / integrate(mps)

		corrs = [v]
		ids2N = [k for k in sampleidx(N) if k > 1]
		idsallN = [1; ids2N]
		for i in ids2N
			pos2 = index(lattice, i, branch=:+)
			m = ADTTerm((pos2,pos1), (zdiag, zdiag))
			mps2 = apply!(m, copy(mps))
			v = integrate(mps2) / integrate(mps)
			push!(corrs, v)
		end
		
		A = kron(op, Ib)
		corrs2 = correlation_2op_1t(H, A, A, ρ, 0:δt:t, reverse = false)
		corrs2 = corrs2[idsallN]
		@test norm(corrs - corrs2) / norm(corrs2) < tol


		# off-diagonal observables

		op1 = [0 0.8; 0 0]
		op2 = [0 0; 0.7*im 0]

		A1 = kron(op1, Ib)
		A2 = kron(op2, Ib)


		c1 = ContourIndex(1, branch=:+)

		ct = ContourOperator(c1, op1 * op2)
		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / integrate(mps)

		corrs = [v]
		c2 = ContourIndex(1, branch=:+)
		for i in ids2N
			c1 = ContourIndex(i, branch=:+)
			ct = ContourOperator([c1, c2], [op1, op2])

			mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
			mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
			mps2 = mult!(mpsK, mpsI, trunc=trunc)
			v = integrate(mps2) / integrate(mps)

			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = false)
		corrs2 = corrs2[idsallN]

		@test norm(corrs - corrs2) / norm(corrs2) < tol

		op1 = [0 0.8*im; 0 0]
		op2 = [0 0; 0.7 0]


		A1 = kron(op1, Ib)
		A2 = kron(op2, Ib)

		c1 = ContourIndex(1, branch=:-)

		ct = ContourOperator(c1, op1 * op2)
		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / integrate(mps)

		corrs = [v]
		for i in ids2N
			c2 = ContourIndex(i, branch=:+)
			ct = ContourOperator([c1, c2], [op1, op2])

			mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
			mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
			mps2 = mult!(mpsK, mpsI, trunc=trunc)
			v = integrate(mps2) / integrate(mps)

			push!(corrs, v)
		end

		corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = true)
		corrs2 = corrs2[idsallN]

		@test norm(corrs - corrs2) / norm(corrs2) < tol

	end

end



@testset "Rabi model: mixed-time" begin

	Ω = 0.5
	Nt = 5
	δt = 0.03
	t = Nt * δt
	Nτ = 10
	δτ = 0.05
	β = Nτ * δτ
	chi = 100
	d = 50
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(Nt = Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed)

	x = [0 1; 1 0]
	xop = Ω .* x
	z = [-1 0; 0 1]
	Is = one(x)
	Ib = one(zeros(d, d))
	model = ImpurityHamiltonian(xop)

	op1 = [0 0.8; 0 0]
	op2 = [0 0; 0.7 0]


	A1 = kron(op1, Ib)
	A2 = kron(op2, Ib)


	bs = AdditiveHyb([z[i,i] for i in 1:size(z,1)])
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI′ = hybriddynamics_naive(lattice, corr, bs, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI′) < tol

	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	mps = mult!(mpsK, mpsI, trunc=trunc)


	H, Hbarebath = rabi_ham(Ω, d=d)


	ρ = exp(-β .* H)

	# off-diagonal observables

	c1 = ContourIndex(1, branch=:+)

	ct = ContourOperator(c1, op1 * op2)
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	c2 = ContourIndex(1, branch=:+)
	ids2Nt = [k for k in sampleidx(Nt) if k > 1]
	idsallNt = [1; ids2Nt]
	for i in ids2Nt
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(H, A1, A2, ρ, 0:δt:t, reverse = false)
	corrs2 = corrs2[idsallNt]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	c1 = ContourIndex(1, branch=:-)

	ct = ContourOperator(c1, op2 * op1)
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	for i in ids2Nt
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op2, op1])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(H, A2, A1, ρ, 0:δt:t, reverse = true)
	corrs2 = corrs2[idsallNt]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	zval = integrate(mps)

	## 虚时格林函数 G(τ) = <op1(τ) op2(0)>：两个算符均插入 τ 分支
	c2 = ContourIndex(1, branch=:τ)
	c1 = ContourIndex(1, branch=:τ)
	ct = ContourOperator(c1, op1 * op2)
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / zval

	corrs = [v]
	ids2Nτ = [k for k in sampleidx(Nτ) if k > 1]
	idsallNτ = [1; ids2Nτ]
	for i in ids2Nτ
		c1 = ContourIndex(i, branch=:τ)
		ct = ContourOperator([c1, c2], [op1, op2])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / zval

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1τ(H, A1, A2, 0:δτ:β, β=β)
	corrs2 = corrs2[idsallNτ]

	@test norm(corrs - corrs2) / norm(corrs2) < tol


	## 虚实混合格林函数 G(τ, t) = <op1(τ) op2(t)>：op1 在 τ 分支，op2 在实时分支
	# ED 参考
	Fed = eigen(Hermitian(H))
	evU, evλ = Fed.vectors, Fed.values
	ρed = evU * Diagonal(exp.(-β .* evλ)) * evU'
	zed = tr(ρed)
	op2ts = [evU * Diagonal(exp.(im .* tj .* evλ)) * evU' * A2 * evU * Diagonal(exp.(-im .* tj .* evλ)) * evU' for tj in 0:δt:t]

	for br in (:+, :-)
		for i in (1, div(Nτ, 2), Nτ)  # sampled τ points (regression)
			τv = (i - 1) * δτ
			op1τ = evU * Diagonal(exp.(τv .* evλ)) * evU' * A1 * evU * Diagonal(exp.(-τv .* evλ)) * evU'
			corrs2 = [tr(op1τ * Bt * ρed) / zed for Bt in op2ts]
			corrs2 = corrs2[1:Nt]

			c1 = ContourIndex(i, branch=:τ)
			c2 = ContourIndex(1, branch=br)
			ct = ContourOperator([c1, c2], [op1, op2])
			mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
			mpsK = boundarycondition!(mpsK, lattice)
			mps2 = mult!(mpsK, mpsI, trunc=trunc)
			v = integrate(mps2) / zval

			corrs = [v]
			for j in 2:Nt
				c2 = ContourIndex(j, branch=br)
				ct = ContourOperator([c1, c2], [op1, op2])

				mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
				mpsK = boundarycondition!(mpsK, lattice)
				mps2 = mult!(mpsK, mpsI, trunc=trunc)
				v = integrate(mps2) / zval

				push!(corrs, v)
			end

			@test norm(corrs - corrs2) / norm(corrs2) < tol
		end
	end

end
println("------------------------------------")
println("|     Dissipative Rabi Model       |")
println("------------------------------------")

@testset "Rabi model: real-time (dissipative)" begin

	Ω = 0.5
	N = 6  # Liouville-space dynamics (d² = 4): keep the chain short
	δt = 0.05
	β = 2
	t = N * δt
	chi = 50
	d = 20
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice = ADTLattice(N = N, δt=δt, contour=:real)

	# x = [0 1; 1 0]
	x = Matrix{ComplexF64}([0 im; -im 0])
	hop = Ω .* x
	z = [-1 0; 0 1]
	Is = one(x)
	Ib = one(zeros(d, d))

	jumpops = [randn(ComplexF64, 2, 2), randn(ComplexF64, 2, 2)]
	model = ImpurityLindbladian(lindbladoperator(hop, jumpops))

	Hbarebath = bosondensityoperator(d=d)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib) + kron(Is, Hbarebath) + kron(z, a' + a)

	jumpops2 = [kron(jump, Ib) for jump in jumpops]
	Lop = lindbladoperator(H, jumpops2)


	bs = AdditiveHyb([z[i,i] for i in 1:size(z,1)])
	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	mpsI = hybriddynamics(lattice, corr, bs, trunc=trunc)
	mpsI′ = hybriddynamics_naive(lattice, corr, bs, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI′) < tol


	ρimp = _rand_dm(2)



	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
	mps = mult!(mpsK, mpsI, trunc=trunc)

	
	ρ = kron(ρimp, exp(-β * Hbarebath)) 

	## diagonal observables
	op = [-0.73 0; 0 0.5]
	zdiag = [op[i,i] for i in 1:size(z, 1)]

	pos1 = index(lattice, 1, branch=:+)
	m = ADTTerm(pos1, zdiag .* zdiag )
	mps2 = apply!(m, copy(mps))
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	ids2N = [k for k in sampleidx(N) if k > 1]
	idsallN = [1; ids2N]
	for i in ids2N
		pos2 = index(lattice, i, branch=:+)
		m = ADTTerm((pos2,pos1), (zdiag, zdiag))
		mps2 = apply!(m, copy(mps))
		v = integrate(mps2) / integrate(mps)
		push!(corrs, v)
	end
	
	A = kron(op, Ib)
	corrs2 = correlation_2op_1t(Lop, A, A, ρ, 0:δt:t, reverse = false)
	corrs2 = corrs2[idsallN]
	@test norm(corrs - corrs2) / norm(corrs2) < tol


	# off-diagonal observables

	op1 = [0 0.8; 0 0]
	op2 = [0 0; 0.7*im 0]

	A1 = kron(op1, Ib)
	A2 = kron(op2, Ib)


	c1 = ContourIndex(1, branch=:+)

	ct = ContourOperator(c1, op1 * op2)
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	c2 = ContourIndex(1, branch=:+)
	for i in ids2N
		c1 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, A1, A2, ρ, 0:δt:t, reverse = false)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

	op1 = [0 0.8*im; 0 0]
	op2 = [0 0; 0.7 0]


	A1 = kron(op1, Ib)
	A2 = kron(op2, Ib)

	c1 = ContourIndex(1, branch=:-)

	ct = ContourOperator(c1, op1 * op2)
	mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
	mps2 = mult!(mpsK, mpsI, trunc=trunc)
	v = integrate(mps2) / integrate(mps)

	corrs = [v]
	for i in ids2N
		c2 = ContourIndex(i, branch=:+)
		ct = ContourOperator([c1, c2], [op1, op2])

		mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
		mps2 = mult!(mpsK, mpsI, trunc=trunc)
		v = integrate(mps2) / integrate(mps)

		push!(corrs, v)
	end

	corrs2 = correlation_2op_1t(Lop, A1, A2, ρ, 0:δt:t, reverse = true)
	corrs2 = corrs2[idsallN]

	@test norm(corrs - corrs2) / norm(corrs2) < tol

end

