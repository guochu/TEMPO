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

# ---- 实时 PT：四支生成元（恢复物理管线验证过的装配约定）----
#
# 与 real.jl 的 pt_ti_mpotensor 调用一一同构（免 Prony 拟合）：
#   · (++): 朴素算符 (op1, op2)——衰减项 op1@晚时配 ηⱼ₊₊ 核、op2@晚时配 ηₖⱼ₊₊ 核，
#     在位 ηⱼ₊₊[1]·(op1·op2) + ηₖⱼ₊₊[1]·(op2·op1)，内容只挂在 + 站点（− 站点直通）。
#   · (−−): 转置算符 (op1ᵀ, op2ᵀ)，同构于 (++)，内容只挂在 − 站点。
#   · (+,−)/(−,+): 跨支项在**时间单元链**（每单元 = fused (b1,b2) 双站点，物理维 d²）
#     上装配——算符为 fused_op(op1, b1)@行 / fused_op(op2, b2)@列，ηⱼₖ/ηₖⱼ 两族，
#     等时项 (ηⱼₖ[1]+ηₖⱼ[1])·(fused_op(op1,b1)·fused_op(op2,b2))；装配后逐单元
#     split 摆位到交错光格（_fit_cross_units_to_lattice，机制同 split_mpotensor）。
# 该约定由 ptmodels 的 freebosons / JC model 实时精确解验证（见 test/ptmodels）。

# 对角支（同支）：全光格 FSM，朴素 d×d 算符，内容挂在 f 支站点
function _exact_diag_generators(lattice::RealPTLattice, η::CorrelationMatrix,
								opA::AbstractMatrix, opB::AbstractMatrix, br::Symbol)
	N, d, L = lattice.N, lattice.d, length(lattice)
	poss(j) = lattice[BosonicContourIndex(j, br)]
	C1, C2 = zeros(ComplexF64, L, L), zeros(ComplexF64, L, L)
	for i in 1:N, j in 1:N
		i > j && (C1[poss(i), poss(j)] = η.ηⱼₖ[i - j + 1])
		i < j && (C2[poss(j), poss(i)] = η.ηₖⱼ[j - i + 1])
	end
	h = fill(zeros(ComplexF64, d, d), L)
	for j in 1:N
		h[poss(j)] = η.ηⱼₖ[1] .* (opA * opB) .+ η.ηₖⱼ[1] .* (opB * opA)
	end
	# fam1 = (opA@晚, opB@早)·ηⱼₖ；fam2 = (opB@晚, opA@早)·ηₖⱼ
	return _exact_branch_ham((opA, opB, C1), (opB, opA, C2), h)
end

# 跨支：时间单元链 FSM（N 单元、时间降序、fused d² 物理腿），fused 行/列算符
function _exact_cross_unit_ham(lattice::RealPTLattice, η::CorrelationMatrix,
							   op1::AbstractMatrix, op2::AbstractMatrix,
							   b1::Symbol, b2::Symbol)
	N, d = lattice.N, lattice.d
	d2 = d * d
	pu(j) = N - j + 1                      # 单元链位（时间降序：晚时在前）
	fA, fB = fused_op(op1, b1), fused_op(op2, b2)
	C1, C2 = zeros(ComplexF64, N, N), zeros(ComplexF64, N, N)
	for i in 1:N, j in 1:N
		i > j && (C1[pu(i), pu(j)] = η.ηⱼₖ[i - j + 1])
		i < j && (C2[pu(j), pu(i)] = η.ηₖⱼ[j - i + 1])
	end
	h = fill(zeros(ComplexF64, d2, d2), N)
	eqt = η.ηⱼₖ[1] + η.ηₖⱼ[1]
	for u in 1:N
		h[u] = eqt .* (fA * fB)
	end
	# fam1 = (fused_op(op1,b1)@晚, fused_op(op2,b2)@早)·ηⱼₖ；
	# fam2 = (fused_op(op2,b2)@晚, fused_op(op1,b1)@早)·ηₖⱼ
	return _exact_branch_ham((fA, fB, C1), (fB, fA, C2), h)
end

function _exact_realtime_generators(lattice::RealPTLattice, corr::RealCorrelationFunction, hyb::GeneralHybStyle)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr)
	op1, op2 = pairop(hyb)
	H⁺⁺ = _exact_diag_generators(lattice, η⁺⁺, op1, op2, :+)
	H⁻⁻ = _exact_diag_generators(lattice, η⁻⁻, transpose(op1), transpose(op2), :-)
	H⁺⁻ = _exact_cross_unit_ham(lattice, η⁺⁻, op1, op2, :+, :-)
	H⁻⁺ = _exact_cross_unit_ham(lattice, η⁻⁺, op1, op2, :-, :+)
	return (H⁺⁺, H⁺⁻, H⁻⁺, H⁻⁻)
end

# 跨支单元链张量 → 交错光格上的 ProcessTensor。与 Prony 路径不同，PairOpTerm
# 的 bulk 张量**逐单元不同**（B 块携带 C[k,s]·b 的站点值），因此不做三张量复制，
# 而是把每个时间单元的 fused d² 张量单独 split：ket 半（含行腿 fused_op(op1,b1)
# 的 ket 扇区）落在 + 站点、bra 半落在 − 站点，单元间键沿链自然衔接。
function _fit_cross_units_to_lattice(lattice::RealPTLattice, ts::AbstractVector,
									 trunc::TruncationScheme=Defaults.alg_trunc_sys())
	L = length(lattice)
	N = lattice.N
	length(ts) == N || throw(DimensionMismatch("need $N unit tensors, got $(length(ts))"))
	data2 = similar(ts, L)
	for j in N:-1:1
		u, v = split_mpotensor(ts[N - j + 1], trunc)
		data2[lattice[BosonicContourIndex(j, :+)]] = u
		data2[lattice[BosonicContourIndex(j, :-)]] = v
	end
	return ProcessTensor(data2)
end

function influenceoperators(lattice::RealPTLattice, corr::RealCorrelationFunction, hyb::GeneralHybStyle, algexpan::ExactConstructor)
	H⁺⁺, H⁺⁻, H⁻⁺, H⁻⁻ = _exact_realtime_generators(lattice, corr, hyb)
	return (ProcessTensor(tompotensors(H⁺⁺)),
	        _fit_cross_units_to_lattice(lattice, tompotensors(H⁺⁻)),
	        _fit_cross_units_to_lattice(lattice, tompotensors(H⁻⁺)),
	        ProcessTensor(tompotensors(H⁻⁻)))
end

function influenceoperatorsteppers(lattice::RealPTLattice, corr::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle,
                                    alg::FirstOrderStepper, algexpan::ExactConstructor)
	H⁺⁺, H⁺⁻, H⁻⁺, H⁻⁻ = _exact_realtime_generators(lattice, corr, hyb)
	return (ProcessTensor(timeevompo(H⁺⁺, dt, alg).data),
	        _fit_cross_units_to_lattice(lattice, timeevompo(H⁺⁻, dt, alg).data),
	        _fit_cross_units_to_lattice(lattice, timeevompo(H⁻⁺, dt, alg).data),
	        ProcessTensor(timeevompo(H⁻⁻, dt, alg).data))
end

function influenceoperatorsteppers(lattice::RealPTLattice, corr::RealCorrelationFunction, dt::Real, hyb::GeneralHybStyle,
                                    alg::ComplexStepper, algexpan::ExactConstructor)
	H⁺⁺, H⁺⁻, H⁻⁺, H⁻⁻ = _exact_realtime_generators(lattice, corr, hyb)
	function pair_diag(H)
		ma, mb = timeevompo(H, dt, alg)
		return ProcessTensor(ma.data), ProcessTensor(mb.data)
	end
	function pair_cross(H)
		ma, mb = timeevompo(H, dt, alg)
		return _fit_cross_units_to_lattice(lattice, ma.data),
		       _fit_cross_units_to_lattice(lattice, mb.data)
	end
	h1a, h1b = pair_diag(H⁺⁺)
	h2a, h2b = pair_cross(H⁺⁻)
	h3a, h3b = pair_cross(H⁻⁺)
	h4a, h4b = pair_diag(H⁻⁻)
	return h1a, h1b, h2a, h2b, h3a, h3b, h4a, h4b
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
