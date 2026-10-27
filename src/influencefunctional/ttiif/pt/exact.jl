# ExactConstructor：跳过指数展开（Prony）的精确影响算符构造
#
# 生成元 Heff 是 all-to-all 两体 + 在位项的算符和：Σ_j h_j + Σ_{a<b} c_{ab} X_a Y_b。
# 构造委托给 FiniteMPSAlgorithms 的 PairOpTerm：同一对（开启算符 a，收尾算符 b）
# 的全部 (i<j) 配对组成一个 PairOpTerm(a, b, C)——每对开启位点共享一条内部通道的
# 有限状态机 MPO（键维 ≤ L+1），构造不含任何近似；同支的第二个通道族（另一对
# 开启/收尾算符）由 MPOHamiltonian 的精确通道直和（FMA `Base.:+`，块对角拼接）
# 合并，在位项挂在第一族的 D 角。用户可在 `ExactConstructor()` 与 Prony 类
# `algexpan` 之间自由切换；与严格参照（如 QuAPI.Dense 的稠密 IF）对齐时建议用
# `ExactConstructor`，避免拟合误差扰动。

"""
    ExactConstructor()

`algexpan` 的精确替代：不做衰减序列的指数展开拟合，而是把影响算符生成元
（all-to-all 两体 + 在位项）按（开启算符，收尾算符）通道族写成 FiniteMPSAlgorithms
的 `PairOpTerm` 有限状态机 MPO（每对开启位点一条共享通道，键维 ≤ L+1），族间用
精确通道直和合并。构造对任意衰减序列（离散谱/连续谱）严格，供精确对齐检验使用。
要求链长 L ≥ 2（`PairOpTerm`/`MPOHamiltonian` 的下限）。
"""
struct ExactConstructor end

# ---- 两通道族 (a, b, C) + 可选每 site 在位项 hlocs 的精确装配 ----
# 在位项只挂第一族的 D 角；第二族无在位项，直和后 D 角保持原值。
function _exact_branch_ham(fam1::Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix},
                           fam2::Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix},
                           hlocs::Union{Nothing, AbstractVector{<:AbstractMatrix}} = nothing)
	H1 = hlocs === nothing ? MPOHamiltonian(PairOpTerm(fam1...)) :
							MPOHamiltonian(PairOpTerm(fam1...), hlocs)
	return H1 + MPOHamiltonian(PairOpTerm(fam2...))
end

# ---- 实时 PT：四支生成元（与修正后的算符约定一致）----

function _exact_realtime_generators(lattice::RealPTLattice, corr::RealCorrelationFunction, hyb::GeneralHybStyle)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr)
	op1, op2 = pairop(hyb)
	N, d, L = lattice.N, lattice.d, length(lattice)
	T1, T2 = transpose(op1), transpose(op2)      # + 支行腿
	A1, A2 = adjoint(op1), adjoint(op2)          # − 支行腿
	poss(b, j) = lattice[BosonicContourIndex(j, b)]
	zc() = zeros(ComplexF64, L, L)

	# (++): 在位 ηⱼ₊₊[1]·(ᵀop2·ᵀop1) + ηₖⱼ₊₊[1]·(ᵀop1·ᵀop2)；
	#       衰减 m1 = (ᵀop1@晚, ᵀop2@早)·ηⱼ₊₊，m2 = (ᵀop2@晚, ᵀop1@早)·ηₖⱼ₊₊
	C⁺¹, C⁺² = zc(), zc()
	for i in 1:N, j in 1:N
		i > j && (C⁺¹[poss(:+, i), poss(:+, j)] = η⁺⁺.ηⱼₖ[i - j + 1])
		i < j && (C⁺²[poss(:+, j), poss(:+, i)] = η⁺⁺.ηₖⱼ[j - i + 1])
	end
	h⁺ = fill(zeros(ComplexF64, d, d), L)
	for j in 1:N
		h⁺[poss(:+, j)] = η⁺⁺.ηⱼₖ[1] .* (T2 * T1) .+ η⁺⁺.ηₖⱼ[1] .* (T1 * T2)
	end
	H⁺⁺ = _exact_branch_ham((T1, T2, C⁺¹), (T2, T1, C⁺²), h⁺)

	# (−−): 在位 ηⱼ₋₋[1]·(op1'·op2') + ηₖⱼ₋₋[1]·(op2'·op1')；衰减 M1↔M2 互换
	C⁻¹, C⁻² = zc(), zc()
	for i in 1:N, j in 1:N
		i > j && (C⁻¹[poss(:-, i), poss(:-, j)] = η⁻⁻.ηⱼₖ[i - j + 1])
		i < j && (C⁻²[poss(:-, j), poss(:-, i)] = η⁻⁻.ηₖⱼ[j - i + 1])
	end
	h⁻ = fill(zeros(ComplexF64, d, d), L)
	for j in 1:N
		h⁻[poss(:-, j)] = η⁻⁻.ηⱼₖ[1] .* (A1 * A2) .+ η⁻⁻.ηₖⱼ[1] .* (A2 * A1)
	end
	H⁻⁻ = _exact_branch_ham((A2, A1, C⁻¹), (A1, A2, C⁻²), h⁻)

	# (+,−) 与 (−,+)：同算符配对（含同 cell 等时），(−,+) 转置读取核。
	# 交错布局下 p⁺ < p⁻ ⇔ i ≥ j：i ≥ j 开启 ᵀop@+ / 收尾 op'@−，i < j 反向。
	C⁺⁻¹, C⁺⁻² = zc(), zc()
	C⁻⁺¹, C⁻⁺² = zc(), zc()
	for i in 1:N, j in 1:N
		c1 = i == j ? η⁺⁻.ηⱼₖ[1] + η⁺⁻.ηₖⱼ[1] : η⁺⁻[i, j]
		c2 = i == j ? η⁻⁺.ηⱼₖ[1] + η⁻⁺.ηₖⱼ[1] : η⁻⁺[j, i]
		p⁺, p⁻ = poss(:+, i), poss(:-, j)
		if p⁺ < p⁻
			C⁺⁻¹[p⁺, p⁻] = c1
			C⁻⁺¹[p⁺, p⁻] = c2
		else
			C⁺⁻²[p⁻, p⁺] = c1
			C⁻⁺²[p⁻, p⁺] = c2
		end
	end
	H⁺⁻ = _exact_branch_ham((T1, A1, C⁺⁻¹), (A1, T1, C⁺⁻²))
	H⁻⁺ = _exact_branch_ham((T2, A2, C⁻⁺¹), (A2, T2, C⁻⁺²))

	return (H⁺⁺, H⁺⁻, H⁻⁺, H⁻⁻)
end

function influenceoperators(lattice::RealPTLattice, corr::RealCorrelationFunction, hyb::GeneralHybStyle, algexpan::ExactConstructor)
	Hs = _exact_realtime_generators(lattice, corr, hyb)
	return map(H -> ProcessTensor(tompotensors(H)), Hs)
end

function influenceoperatorsteppers(lattice::RealPTLattice, corr::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle,
                                    alg::FirstOrderStepper, algexpan::ExactConstructor)
	Hs = _exact_realtime_generators(lattice, corr, hyb)
	return map(H -> ProcessTensor(timeevompo(H, dt, alg).data), Hs)
end

function influenceoperatorsteppers(lattice::RealPTLattice, corr::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle,
                                    alg::ComplexStepper, algexpan::ExactConstructor)
	Hs = _exact_realtime_generators(lattice, corr, hyb)
	return ntuple(8) do k
		H = Hs[(k + 1) ÷ 2]
		m = timeevompo(H, dt, alg)
		ProcessTensor((isodd(k) ? first(m) : last(m)).data)
	end
end

function influenceoperatorstepper(lattice::RealPTLattice, corr::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle,
                                   alg::FirstOrderStepper, algmult::MPSAlgorithm, algexpan::ExactConstructor)
	h1, h2, h3, h4 = influenceoperatorsteppers(lattice, corr, dt, hyb, alg, algexpan)
	mps = mult(h2, h1, algmult)
	mps = mult(h3, mps, algmult)
	mps = mult(h4, mps, algmult)
	return mps
end

function influenceoperatorstepper(lattice::RealPTLattice, corr::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle,
                                   alg::ComplexStepper, algmult::MPSAlgorithm, algexpan::ExactConstructor)
	h1a, h1b, h2a, h2b, h3a, h3b, h4a, h4b = influenceoperatorsteppers(lattice, corr, dt, hyb, alg, algexpan)
	mps = mult(h1b, h1a, algmult)
	mps = mult(h2a, mps, algmult)
	mps = mult(h2b, mps, algmult)
	mps = mult(h3a, mps, algmult)
	mps = mult(h3b, mps, algmult)
	mps = mult(h4a, mps, algmult)
	mps = mult(h4b, mps, algmult)
	return mps
end

# ---- 虚时 PT：单支生成元（plain 算符，与 imag.jl 约定一致）----

function _exact_imag_generators(lattice::ImagPTLattice, corr::ImagCorrelationFunction, hyb::GeneralHybStyle)
	op1, op2 = pairop(hyb)
	N, d, L = lattice.N, lattice.d, length(lattice)
	η = corr.data
	poss(j) = lattice[BosonicContourIndex(j, :τ)]
	C1, C2 = zeros(ComplexF64, L, L), zeros(ComplexF64, L, L)
	for i in 1:N, j in 1:N
		i > j && (C1[poss(i), poss(j)] = η.ηⱼₖ[i - j + 1])
		i < j && (C2[poss(j), poss(i)] = η.ηₖⱼ[j - i + 1])
	end
	h = fill(zeros(ComplexF64, d, d), L)
	for j in 1:N
		h[poss(j)] = η.ηⱼₖ[1] .* (op1 * op2) .+ η.ηₖⱼ[1] .* (op2 * op1)
	end
	return (_exact_branch_ham((op1, op2, C1), (op2, op1, C2), h),)
end

function influenceoperators(lattice::ImagPTLattice, corr::ImagCorrelationFunction, hyb::GeneralHybStyle, algexpan::ExactConstructor)
	Hs = _exact_imag_generators(lattice, corr, hyb)
	return map(H -> ProcessTensor(tompotensors(H)), Hs)
end

function influenceoperatorsteppers(lattice::ImagPTLattice, corr::ImagCorrelationFunction, dt::Real, hyb::GeneralHybStyle,
                                    alg::FirstOrderStepper, algexpan::ExactConstructor)
	Hs = _exact_imag_generators(lattice, corr, hyb)
	return map(H -> ProcessTensor(timeevompo(H, dt, alg).data), Hs)
end

function influenceoperatorsteppers(lattice::ImagPTLattice, corr::ImagCorrelationFunction, dt::Real, hyb::GeneralHybStyle,
                                    alg::ComplexStepper, algexpan::ExactConstructor)
	(H,) = _exact_imag_generators(lattice, corr, hyb)
	ma, mb = timeevompo(H, dt, alg)
	return (ProcessTensor(ma.data), ProcessTensor(mb.data))
end

function influenceoperatorstepper(lattice::ImagPTLattice, corr::ImagCorrelationFunction, dt::Real, hyb::GeneralHybStyle,
                                   alg::FirstOrderStepper, algmult::MPSAlgorithm, algexpan::ExactConstructor)
	(mps,) = influenceoperatorsteppers(lattice, corr, dt, hyb, alg, algexpan)
	return mps
end

function influenceoperatorstepper(lattice::ImagPTLattice, corr::ImagCorrelationFunction, dt::Real, hyb::GeneralHybStyle,
                                   alg::ComplexStepper, algmult::MPSAlgorithm, algexpan::ExactConstructor)
	mps1, mps2 = influenceoperatorsteppers(lattice, corr, dt, hyb, alg, algexpan)
	return mult(mps2, mps1, algmult)
end
