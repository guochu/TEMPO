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
# ALS 收敛后 bra 链本身就是最终输出态，规范与 Schmidt 谱不再依赖环境：直接用
# FiniteMPSAlgorithms 的 `canonicalize!`——先 QR 左扫（精确、不截断、不改变
# 状态），再以 `trunc` 做 SVD 右扫（每步把 u·s 折回左侧并把谱写入该键的
# Schmidt 值）。两步扫之后 bra 右正交，`z.s` 是（截断后）状态在各键的精确
# Schmidt 值；这一保证不依赖 ALS 是否到达不动点。
# ---------------------------------------------------------------------------
function _finalize!(m::Union{MultCache,HadamardCache}, alg::FiniteMPSAlgorithms.DMRG1, trunc::TruncationScheme=NoTruncation())
	canonicalize!(m.bra; alg=Orthogonalize(SVD(), trunc))
	return m
end
