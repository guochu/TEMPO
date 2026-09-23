# ---------------------------------------------------------------------------
# FiniteMPSAlgorithms 适配层
#
# ADT / ProcessTensor 的存储 payload 直接内嵌 FiniteMPSAlgorithms 的
# CanonicalMPS / CanonicalMPO（字段 `.data`，见 adt/def.jl 与 pt/def.jl）：
# FiniteMPSAlgorithms 的算法直接作用在 payload 上并就地改写，TEMPO 侧无需
# 任何转换。此外提供 TEMPO `DMRG1` 配置到 FMA 迭代参数的翻译。
# ---------------------------------------------------------------------------

# TEMPO 的 `DMRG1` 携带截断方案（`trunc.D` 是键维上限，并用于 initguess 与
# finalize 截断）以及 initguess/callback 字段；FiniteMPSAlgorithms 的 `DMRG1`
# 只有纯迭代参数（maxiter/tol/D/verbosity）。此翻译器把 TEMPO 配置映射到
# ALS 引擎的迭代参数。
_fmadmrg1(alg::DMRG1) = FiniteMPSAlgorithms.DMRG1(maxiter=alg.maxiter, tol=alg.tol, D=alg.trunc.D, verbosity=alg.verbosity)
