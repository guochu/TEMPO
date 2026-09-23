# ---------------------------------------------------------------------------
# FiniteMPSAlgorithms 适配层
#
# ADT / ProcessTensor 的存储 payload 直接内嵌 FiniteMPSAlgorithms 的
# CanonicalMPS / CanonicalMPO（字段 `.parent`，见 adt/def.jl 与 pt/def.jl）：
# FiniteMPSAlgorithms 的算法直接作用在 payload 上并就地改写，TEMPO 侧无需
# 任何转换。此外提供 TEMPO `DMRG1` 配置到 FMA 迭代参数的翻译，以及
# ADT/PT 的 `mult` 共用的 finalize sweep。
# ---------------------------------------------------------------------------

# TEMPO 的 `DMRG1` 携带截断方案（`trunc.D` 是键维上限，并用于 initguess 与
# finalize 截断）以及 initguess/callback 字段；FiniteMPSAlgorithms 的 `DMRG1`
# 只有纯迭代参数（maxiter/tol/D/verbosity）。此翻译器把 TEMPO 配置映射到
# ALS 引擎的迭代参数。
_fmadmrg1(alg::DMRG1) = FiniteMPSAlgorithms.DMRG1(maxiter=alg.maxiter, tol=alg.tol, D=alg.trunc.D, verbosity=alg.verbosity)

# ---------------------------------------------------------------------------
# finalize sweep（ADT/PT 的 `mult`（DMRG1 路线）共用）
#
# 与旧版 `rightsweep_final!` 一致：QR 左扫重建左环境后，从右到左以 `trunc`
# 做（可截断的）SVD 重新规范化，并把归一化键谱写入 Schmidt 值。
# PT 侧以 `NoTruncation()` 调用（保持旧行为：只重新规范化 + 写谱）。
# ---------------------------------------------------------------------------
function _finalize!(m::MultCache, alg::FiniteMPSAlgorithms.DMRG1, trunc::TruncationScheme=NoTruncation())
	L = length(m.bra)
	leftsweep!(m, alg)
	for site in L:-1:2
		mpsj = _reduce_site(m.ket[site], m.H[site], m.hstorage[site], m.hstorage[site+1])
		u, s, v = tsvd!(mpsj, (1,), (2, 3, 4), trunc=trunc)
		m.bra[site] = v
		if site == 2
			m.bra[1] = _contract_last(m.bra[1], u * Diagonal(s))
		end
		m.bra.s[site] = normalize!(s)
		m.hstorage[site] = _env_updateright(m.hstorage[site+1], m.bra[site], m.H[site], m.ket[site])
	end
	return m
end

# HadamardCache 版本（ADT×ADT 的 mult）：rank-3 site 张量、物理指标共享
function _finalize!(m::HadamardCache, alg::FiniteMPSAlgorithms.DMRG1, trunc::TruncationScheme=NoTruncation())
	L = length(m.bra)
	leftsweep!(m, alg)
	for site in L:-1:2
		mpsj = _reduce_hadamard_site(m.ketx[site], m.kety[site], m.hstorage[site], m.hstorage[site+1])
		u, s, v = tsvd!(mpsj, (1,), (2, 3), trunc=trunc)
		m.bra[site] = permute(v, (1, 2), (3,))
		if site == 2
			m.bra[1] = _contract_last(m.bra[1], u * Diagonal(s))
		end
		m.bra.s[site] = normalize!(s)
		m.hstorage[site] = _updateright(m.hstorage[site+1], m.bra[site], m.ketx[site], m.kety[site])
	end
	return m
end
