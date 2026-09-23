# ---------------------------------------------------------------------------
# 算法配置（后端：FiniteMPSAlgorithms）
#
# `Orthogonalize`、`SVDCompression` 与张量层原语（TruncationScheme、tsvd!、
# leftorth!、QR/SVD…）直接采用 FiniteMPSAlgorithms 的类型与实现；
# `MPSAlgorithm` 也来自 FiniteMPSAlgorithms。
#
# `DMRG1` 保留 TEMPO 的公共接口：它携带截断方案 `trunc`（`trunc.D` 既是
# 初始猜测的键维上限，也用于 ALS 收敛后的 finalize 截断）、初始猜测方式
# `initguess` 与 `callback` 字段。FiniteMPSAlgorithms 的 `DMRG1` 只有纯迭代
# 参数（maxiter/tol/D/verbosity，ALS 本身不截断）；两者的翻译与 TEMPO 特有
# 的 finalize（带截断的末次 sweep，见 adt/mult 与 pt/mult）在 wrapper 层完成。
# ---------------------------------------------------------------------------

abstract type DMRGAlgorithm <: MPSAlgorithm end

# truncation schemes carrying an explicit maximum bond dimension `D`; `DMRG1`
# requires one of these, since `D` seeds the initial guess of the sweeps
const TruncationWithD = Union{TruncateDim, TruncateDimCutoff}

const AllowedInitGuesses = (:svd, :pre, :rand)

"""
	DMRG1 <: DMRGAlgorithm

Configuration of an MPS/MPO product compression algorithm based on DMRG iterative sweeping.

# Fields
- `trunc::TruncationWithD`: truncation scheme carrying a maximum bond dimension `D` (`truncdim(D)` or `truncdimcutoff(D, ϵ)`), used both to compress the result and to seed the initial guess
- `maxiter::Int`: maximum number of iterations
- `tol::Float64`: convergence tolerance
- `initguess::Symbol`: initial guess, one of `:svd`, `:pre`, `:rand`
- `verbosity::Int`: verbosity level
- `callback::Function`: callback function

Main constructor: `DMRG1(trunc; maxiter=5, tol=1e-12, initguess=:svd, verbosity=0, callback=Returns(nothing))`.
"""
struct DMRG1{T<:TruncationWithD} <: DMRGAlgorithm
	trunc::T
	maxiter::Int
	tol::Float64
	initguess::Symbol
	verbosity::Int
	callback::Function
end
"""
	DMRG1(trunc::TruncationWithD; maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0, callback::Function=Returns(nothing))

Construct a `DMRG1` algorithm configuration.

# Arguments
- `trunc::TruncationWithD`: truncation scheme carrying a maximum bond dimension `D` (constructed with `truncdimcutoff(D, ϵ)` or `truncdim(D)`); schemes without a dimension cap (`truncrelerr`, `NoTruncation`) are not allowed
- `maxiter::Int`: maximum number of iterations
- `tol::Float64`: convergence tolerance
- `initguess::Symbol`: initial guess, must be one of `:svd`, `:pre`, `:rand`, otherwise an `ArgumentError` is thrown
- `verbosity::Int`: verbosity level
- `callback::Function`: callback function
"""
function DMRG1(trunc::TruncationWithD; maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0, callback::Function=Returns(nothing))
	(initguess in AllowedInitGuesses) || throw(ArgumentError("initguess must be one of $(AllowedInitGuesses)"))
	return DMRG1(trunc, maxiter, tol, initguess, verbosity, callback)
end
"""
	DMRG1(; trunc::TruncationWithD=DefaultITruncation, kwargs...)

Construct a `DMRG1` from keyword arguments, with default truncation scheme `DefaultITruncation`.
"""
DMRG1(; trunc::TruncationWithD=DefaultITruncation, kwargs...) = DMRG1(trunc; kwargs...)
Base.similar(x::DMRG1; trunc::TruncationWithD=x.trunc, maxiter::Int=x.maxiter, tol::Float64=x.tol, initguess::Symbol=x.initguess, verbosity::Int=x.verbosity, callback=x.callback) = DMRG1(
			trunc=trunc, maxiter=maxiter, tol=tol, initguess=initguess, verbosity=verbosity, callback=callback)

# ---------------------------------------------------------------------------
# `SVDCompression`：FiniteMPSAlgorithms 的类型；这里补充 TEMPO 风格的
# positional 构造器（`SVDCompression(trunc; verbosity)`），保持与旧版 TEMPO
# 完全兼容。关键字构造（`SVDCompression(; trunc, verbosity)`）沿用
# FiniteMPSAlgorithms 的 `@kwdef` 定义。
# ---------------------------------------------------------------------------
"""
    SVDCompression(trunc::TruncationScheme; verbosity=0)
    SVDCompression(; trunc=truncdimcutoff(D=Defaults.D, ϵ=Defaults.tol, add_back=0), verbosity=0)

Parameters for an SVD-based DMRG compression algorithm: singular values are truncated according to the truncation scheme `trunc` (any `TruncationScheme`),
while `verbosity` controls the verbosity of the output. The scheme is accessible through the `trunc` field.
"""
SVDCompression(trunc::TruncationScheme; verbosity::Int=0) = SVDCompression{typeof(trunc)}(trunc, verbosity)
Base.similar(x::SVDCompression; trunc::TruncationScheme=x.trunc, verbosity::Int=x.verbosity) = SVDCompression(trunc; verbosity=verbosity)
