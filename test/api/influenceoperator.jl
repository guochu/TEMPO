# influenceoperators / influenceoperatorsteppers / influenceoperatorstepper vs manual
# accumulation of the influence operator, parametrized over ADT (ADTTerm) and PT (FockTermS)

using QuAPI
using QuAPI: bosonic_Δt, BosonicContourIndex
using QuAPI.Dense: BosonicIF, BosonicHybridizationMatrix, BosonicSplitDelta
using FiniteMPSAlgorithms: NoTruncation

@testset "InfluenceOperator: ADT lattice" begin
	# --- imaginary time ---
	β = 1
	δτ = 0.2
	tol = 1.0e-6
	N = round(Int, β/δτ)
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=50, ϵ=1.0e-6, add_back=0)
	algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
	d = 2

	lattice = ADTLattice(N=N, δτ=δτ, d=d, contour=:imag)
	corr = correlationfunction(bath, lattice)
	hyb = AdditiveHyb(randn(Float64, d))

	z = hyb.op
	z2 = z .* z
	zz = reshape(kron(z, z), d, d)

	mpo1 = only(influenceoperators(lattice, corr, hyb, algexpan))

	orth = Orthogonalize(SVD(), trunc)
	mpo2 = nothing
	for i in 1:lattice.N, j in 1:lattice.N
		ind1, ind2 = BosonicContourIndex(i+1), BosonicContourIndex(j+1)
		coef = corr[BosonicContourIndex(i, :τ), BosonicContourIndex(j, :τ)]
		if coef != 0
			if ind1 == ind2
				t = ADTTerm(lattice[ind1], coef .* z2)
			else
				t = ADTTerm((lattice[ind1], lattice[ind2]), coef .* zz)
			end
			if isnothing(mpo2)
				mpo2 = apply!(t, vacuumstate(Float64, lattice))
			else
				mpo2 += apply!(t, vacuumstate(Float64, lattice))
			end
			canonicalize!(mpo2, alg=orth)
		end
	end
	@test distance(mpo1, mpo2) / norm(mpo2) < tol

	# one first-order step: I + dt·K
	dt = 0.01
	mpo1 = dt * mpo1 + vacuumstate(Float64, lattice)
	mps0, = influenceoperatorsteppers(lattice, corr, dt, hyb, WII(), algexpan)
	@test distance(mpo1, mps0) / norm(mps0) < dt

	# single steppers for all four stepper types and two multiplication algorithms
	for algmult in (SVDCompression(trunc=truncdimcutoff(D=50, ϵ=1.0e-12, add_back=0)), DMRG1(trunc=truncdimcutoff(D=50, ϵ=1.0e-6)))
		mps1 = influenceoperatorstepper(lattice, corr, dt, hyb, WII(), algmult, algexpan)
		_n = norm(mps1)
		mps2 = influenceoperatorstepper(lattice, corr, dt, hyb, WI(), algmult, algexpan)
		mps3 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WI()), algmult, algexpan)
		mps4 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WII()), algmult, algexpan)
		@test distance(mps1, mps0) / _n < dt
		@test distance(mps2, mps0) / _n < dt
		@test distance(mps3, mps0) / _n < dt
		@test distance(mps4, mps0) / _n < dt
	end

	# --- real time: one partial influence operator per branch combination ---
	β = 1.0
	δt = 0.1
	N = 4
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-8, add_back=0)
	algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
	tol = 1.0e-5

	lattice = ADTLattice(N=N, δt=δt, contour=:real)
	corr = correlationfunction(bath, lattice)
	hyb = AdditiveHyb(randn(Float64, d))
	z = hyb.op
	z2 = z .* z
	zz = reshape(kron(z, z), d, d)

	function manual_if(b1, b2)
		mps = nothing
		for i in 1:lattice.N, j in 1:lattice.N
			ind1, ind2 = BosonicContourIndex(i, b1), BosonicContourIndex(j, b2)
			coef = corr[BosonicContourIndex(i, b1), BosonicContourIndex(j, b2)]
			if coef != 0
				if ind1 == ind2
					t = ADTTerm(lattice[ind1], coef .* z2)
				else
					t = ADTTerm((lattice[ind1], lattice[ind2]), coef .* zz)
				end
				if isnothing(mps)
					mps = apply!(t, vacuumstate(ComplexF64, lattice))
				else
					mps += apply!(t, vacuumstate(ComplexF64, lattice))
				end
			end
		end
		return mps
	end

	ifs = influenceoperators(lattice, corr, hyb, algexpan)
	@test length(ifs) == 4
	for (b1, b2) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
		mps_ref = manual_if(b1, b2)
		mps = (b1 == :+ && b2 == :+) ? ifs[1] : (b1 == :+ ? ifs[2] : (b2 == :+ ? ifs[3] : ifs[4]))
		@test distance(mps, mps_ref) / norm(mps_ref) < tol
	end

	# steppers: one first-order step for each branch combination
	dt = 0.01
	mps_ref = [dt * m + vacuumstate(ComplexF64, lattice) for m in ifs]
	mps_stp = influenceoperatorsteppers(lattice, corr, dt, hyb, WII(), algexpan)
	for (m1, m2) in zip(mps_ref, mps_stp)
		@test distance(m1, m2) / norm(m2) < dt
	end

	mps0 = mult!(mps_stp[1], mps_stp[2], trunc=trunc)
	mps0 = mult!(mps0, mps_stp[3], trunc=trunc)
	mps0 = mult!(mps0, mps_stp[4], trunc=trunc)
	for algmult in (SVDCompression(trunc=truncdimcutoff(D=50, ϵ=1.0e-12, add_back=0)), DMRG1(trunc=truncdimcutoff(D=50, ϵ=1.0e-6)))
		mps1 = influenceoperatorstepper(lattice, corr, dt, hyb, WII(), algmult, algexpan)
		_n = norm(mps1)
		mps2 = influenceoperatorstepper(lattice, corr, dt, hyb, WI(), algmult, algexpan)
		mps3 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WI()), algmult, algexpan)
		mps4 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WII()), algmult, algexpan)
		@test distance(mps1, mps0) / _n < dt
		@test distance(mps2, mps0) / _n < dt
		@test distance(mps3, mps0) / _n < dt
		@test distance(mps4, mps0) / _n < dt
	end
end

@testset "InfluenceOperator: PT lattice" begin
	__ti_localop(η, op1, op2) = η.ηⱼₖ[1] .* (op1 * op2) + η.ηₖⱼ[1] .* (op2 * op1)

	# --- imaginary time ---
	β = 1
	δτ = 0.2
	tol = 1.0e-6
	N = round(Int, β/δτ)
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=50, ϵ=1.0e-6, add_back=0)
	algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
	d = 2

	lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)
	corr = correlationfunction(bath, lattice)

	for hyb in (NonAdditiveHyb(_rand_ham(d)), NonDiagonalHyb(randn(ComplexF64, d, d)))
		op1, op2 = pairop(hyb)

		mpo1 = only(influenceoperators(lattice, corr, hyb, algexpan))

		mpo2 = nothing
		for i in 1:lattice.N, j in 1:lattice.N
			ind1, ind2 = BosonicContourIndex(i), BosonicContourIndex(j)
			if ind1 == ind2
				m = __ti_localop(corr.data, op1, op2)
				t = FockTermS(lattice[ind1], m)
			else
				coef = corr[BosonicContourIndex(i, :τ), BosonicContourIndex(j, :τ)]
				m = coef .* kron(op2, op1)
				t = FockTermS((lattice[ind1], lattice[ind2]), reshape(m, (d, d, d, d)))
			end
			if isnothing(mpo2)
				mpo2 = apply!(t, vacuumstate(scalartype(hyb), lattice))
			else
				mpo2 += apply!(t, vacuumstate(scalartype(hyb), lattice))
			end
		end
		canonicalize!(mpo2, alg=Orthogonalize(trunc=trunc))
		@test distance(mpo1, mpo2) / norm(mpo2) < tol

		dt = 0.01
		mpo1 = dt * mpo1 + vacuumstate(Float64, lattice)
		mps0, = influenceoperatorsteppers(lattice, corr, dt, hyb, WII(), algexpan)
		@test distance(mpo1, mps0) / norm(mps0) < dt

		for algmult in (SVDCompression(trunc=truncdimcutoff(D=50, ϵ=1.0e-12, add_back=0)), DMRG1(trunc=truncdimcutoff(D=50, ϵ=1.0e-6)))
			mps1 = influenceoperatorstepper(lattice, corr, dt, hyb, WII(), algmult, algexpan)
			_n = norm(mps1)
			mps2 = influenceoperatorstepper(lattice, corr, dt, hyb, WI(), algmult, algexpan)
			mps3 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WI()), algmult, algexpan)
			mps4 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WII()), algmult, algexpan)
			@test distance(mps1, mps0) / _n < dt
			@test distance(mps2, mps0) / _n < dt
			@test distance(mps3, mps0) / _n < dt
			@test distance(mps4, mps0) / _n < dt
		end
	end

	# --- real time ---
	# 手工参照与修正后的构造约定一致（对严格 IF 验证，见 quapi_crosscheck）：
	#  · (++)：衰减项行腿为 (opᵀ, opᵀ)，i>j 用 ηⱼₖ 核配 (ᵀop1@i, ᵀop2@j)，i<j 用 ηₖⱼ 核配
	#    (ᵀop2@i, ᵀop1@j)；等时矩阵积用原算符 (op1·op2)/(op2·op1)。
	#  · (−−)：衰减项行腿为 (op', op') 且 M1↔M2 互换；等时积用 (op1ᵀ·op2ᵀ)/(op2ᵀ·op1ᵀ)。
	#  · (+,−)/(−,+)：同算符配对——(+,−) 配 η⁺⁻ 完整块 (ᵀop1@+, op1'@−)，
	#    (−,+) 配 η⁻⁺ 转置读取 (op2'@−, ᵀop2@+)；等时为求和核配同形式双格点积。
	__get_contour_op(lattice, ind1, ind2, op1, op2, corr) = begin
		d = lattice.d
		b1, b2 = branch(ind1), branch(ind2)
		η = branch(corr, b1, b2)
		pos1, pos2 = lattice[ind1], lattice[ind2]
		if (b1, b2) == (:+, :+)
			if ind1 == ind2
				m = __ti_localop(η, op1, op2)
				t = FockTermS(pos1, m)
			else
				coef = η[ind1.j, ind2.j]
				m = ind1.j > ind2.j ? coef .* kron(transpose(op2), transpose(op1)) :
									  coef .* kron(transpose(op1), transpose(op2))
				t = FockTermS((pos1, pos2), reshape(m, (d, d, d, d)))
			end
		elseif (b1, b2) == (:-, :-)
			if ind1 == ind2
				m = __ti_localop(η, transpose(op1), transpose(op2))
				t = FockTermS(pos1, m)
			else
				coef = η[ind1.j, ind2.j]
				m = ind1.j > ind2.j ? coef .* kron(adjoint(op1), adjoint(op2)) :
									  coef .* kron(adjoint(op2), adjoint(op1))
				t = FockTermS((pos1, pos2), reshape(m, (d, d, d, d)))
			end
		elseif (b1, b2) == (:+, :-)
			coef = ind1 == ind2 ? η.ηⱼₖ[1] + η.ηₖⱼ[1] : η[ind1.j, ind2.j]
			m = coef .* kron(adjoint(op1), transpose(op1))
			t = FockTermS((pos1, pos2), reshape(m, (d, d, d, d)))
		else
			coef = ind1 == ind2 ? η.ηⱼₖ[1] + η.ηₖⱼ[1] : η[ind2.j, ind1.j]
			m = coef .* kron(transpose(op2), adjoint(op2))
			t = FockTermS((pos1, pos2), reshape(m, (d, d, d, d)))
		end
		return t
	end

	β = 1.0
	δt = 0.1
	N = 4
	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β, μ=0)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-8, add_back=0)
	algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
	tol = 1.0e-5

	lattice = PTLattice(N=N, δt=δt, contour=:real)
	corr = correlationfunction(bath, lattice)

	for hyb in (NonAdditiveHyb(_rand_ham(d)), NonDiagonalHyb(randn(ComplexF64, d, d)))
		op1, op2 = pairop(hyb)

		function manual_if(b1, b2)
			mps = nothing
			for i in 1:lattice.N, j in 1:lattice.N
				ind1, ind2 = BosonicContourIndex(i, b1), BosonicContourIndex(j, b2)
				t = __get_contour_op(lattice, ind1, ind2, op1, op2, corr)
				if isnothing(mps)
					mps = apply!(t, vacuumstate(scalartype(hyb), lattice))
				else
					mps += apply!(t, vacuumstate(scalartype(hyb), lattice))
				end
			end
			return mps
		end

		# 非对角耦合的行腿/配对约定用 QuAPI.Dense 的（已对严格 IF 验证的）PT IF 直接对齐检验，
	# 单模 bath + 未截断：若结构正确，两者应在 Prony 拟合误差内一致。
	if hyb isa NonDiagonalHyb
		ω, β0, δt0, N0, Gc = 1.3, 1.0, 0.1, 3, 0.05
		spec0 = DiscreteSpectrum([ω], [Gc^2])
		bath0 = bosonicbath(spec0, β=β0)
		lat0 = PTLattice(N=N0, δt=δt0, contour=:real, d=d)
		corr0 = correlationfunction(bath0, lat0)
		IF0 = influenceoperatorstepper(lat0, corr0, 1.0, hyb, WII(),
			SVDCompression(trunc=NoTruncation()))
		# TEMPO 张量（时间降序 (+,−)/时刻）→ QuAPI 交错布局（时间升序 (+,−)/时刻）
		M0 = TEMPO.todense(IF0.parent; order=:lsb)
		L0 = length(IF0)
		T0 = permutedims(reshape(M0, ntuple(_->d, L0)..., ntuple(_->d, L0)...),
			vcat([[p, L0+p] for p in 1:L0]...))
		perir = vcat([[2t-1, 2t] for t in N0:-1:1]...)
		AT0 = permutedims(T0, vcat([[2p-1, 2p] for p in perir]...))
		corrq0 = bosonic_Δt(spec0; β=β0, N=N0, t=N0*δt0)
		hq0 = BosonicHybridizationMatrix(corrq0; split=true)
		keep0 = trues(2*(N0+1)); keep0[1] = keep0[2] = false
		idx0 = findall(keep0)
		basisq0 = BosonicContourIndex[]
		for k in 1:N0, b in (:+, :-); push!(basisq0, BosonicContourIndex(k+1, b)); end
		h2q0 = BosonicHybridizationMatrix(basisq0,
			QuAPI.Dense.BosonicSplitDelta(hq0.Δ.M1[idx0, idx0], hq0.Δ.M2[idx0, idx0]))
		Tq0 = reshape(BosonicIF(h2q0, op1).data, ntuple(_->d, 4N0))
		@test maximum(abs.(AT0 .- Tq0)) < 1e-6
	end

	ifs = influenceoperators(lattice, corr, hyb, algexpan)
	@test length(ifs) == 4

		dt = 0.01
		mps_ref = [dt * m + vacuumstate(ComplexF64, lattice) for m in ifs]
		mps_stp = influenceoperatorsteppers(lattice, corr, dt, hyb, WII(), algexpan)
		for (m1, m2) in zip(mps_ref, mps_stp)
			@test distance(m1, m2) / norm(m2) < dt
		end

		mps0 = mult!(mps_stp[1], mps_stp[2], trunc=trunc)
		mps0 = mult!(mps0, mps_stp[3], trunc=trunc)
		mps0 = mult!(mps0, mps_stp[4], trunc=trunc)
		for algmult in (SVDCompression(trunc=truncdimcutoff(D=50, ϵ=1.0e-12, add_back=0)), DMRG1(trunc=truncdimcutoff(D=50, ϵ=1.0e-6)))
			mps1 = influenceoperatorstepper(lattice, corr, dt, hyb, WII(), algmult, algexpan)
			_n = norm(mps1)
			mps2 = influenceoperatorstepper(lattice, corr, dt, hyb, WI(), algmult, algexpan)
			mps3 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WI()), algmult, algexpan)
			mps4 = influenceoperatorstepper(lattice, corr, dt, hyb, ComplexStepper(WII()), algmult, algexpan)
			@test distance(mps1, mps0) / _n < dt
			@test distance(mps2, mps0) / _n < dt
			@test distance(mps3, mps0) / _n < dt
			@test distance(mps4, mps0) / _n < dt
		end
	end
end

@testset "InfluenceOperator: ExactConstructor vs Prony" begin
	# 单模离散谱 bath：各支 lag 序列是有限个纯指数的和（正温下 k± 两个模），
	# DeterminedProny 在无噪数据上按阶自动收敛到机器精度残差——这正是
	# "Prony 分解误差很小"的情形。此时 ExactConstructor（无拟合）与 Prony
	# 后端的 influenceoperators 必须张量级一致（稠密链矩阵相对偏差 < 1e-9）。
	d = 2
	β = 1.0
	spec = DiscreteSpectrum([1.3], [0.05])
	bath = bosonicbath(spec, β=β)
	exact = ExactConstructor()
	prony = DeterminedProny(n=4, tol=1.0e-13, verbosity=0)
	tol = 1.0e-9

	_dense(x::ProcessTensor) = TEMPO.todense(x.parent; order=:lsb)

	# 先显式确认 Prony 残差确实很小（对实时 (+,+) 支 ηⱼₖ lag 序列）
	seq = branch(correlationfunction(bath, PTLattice(N=5, δt=0.1, contour=:real, d=d)), :+, :+).ηⱼₖ[2:end]
	xs, λs = exponential_expansion(seq, prony)
	@test expansion_error(seq, xs, λs) < 1.0e-12

	# --- real time：四支逐支对比 ---
	N, δt = 5, 0.1
	lat = PTLattice(N=N, δt=δt, contour=:real, d=d)
	corr = correlationfunction(bath, lat)
	for hyb in (NonAdditiveHyb(_rand_ham(d)), NonDiagonalHyb(randn(ComplexF64, d, d)))
		Iex = influenceoperators(lat, corr, hyb, exact)
		Ipr = influenceoperators(lat, corr, hyb, prony)
		@test length(Iex) == length(Ipr) == 4
		for (m1, m2) in zip(Iex, Ipr)
			D1, D2 = _dense(m1), _dense(m2)
			@test norm(D1 - D2) / norm(D2) < tol
		end
	end

	# --- imaginary time ---
	Nτ, δτ = 5, 0.2
	latτ = PTLattice(N=Nτ, δτ=δτ, contour=:imag, d=d)
	corrτ = correlationfunction(bath, latτ)
	for hyb in (NonAdditiveHyb(_rand_ham(d)), NonDiagonalHyb(randn(ComplexF64, d, d)))
		m1 = only(influenceoperators(latτ, corrτ, hyb, exact))
		m2 = only(influenceoperators(latτ, corrτ, hyb, prony))
		D1, D2 = _dense(m1), _dense(m2)
		@test norm(D1 - D2) / norm(D2) < tol
	end

	# --- steppers 的一致性（一阶步进，实时）；WII 级数演化舍入更大，放宽到 1e-8 ---
	lat3 = PTLattice(N=3, δt=0.1, contour=:real, d=d)
	corr3 = correlationfunction(bath, lat3)
	hyb3 = NonDiagonalHyb(randn(ComplexF64, d, d))
	dt = 0.01
	Sex = influenceoperatorsteppers(lat3, corr3, dt, hyb3, WII(), exact)
	Spr = influenceoperatorsteppers(lat3, corr3, dt, hyb3, WII(), prony)
	for (m1, m2) in zip(Sex, Spr)
		D1, D2 = _dense(m1), _dense(m2)
		@test norm(D1 - D2) / norm(D2) < 1.0e-8
	end
end
