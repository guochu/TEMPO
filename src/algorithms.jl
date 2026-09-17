abstract type MPSAlgorithm end
abstract type DMRGAlgorithm <: MPSAlgorithm end



"""
    SVDCompression(trunc::TruncationScheme; verbosity=0)
    SVDCompression(; trunc=truncdimcutoff(D=Defaults.D, ϵ=Defaults.tol, add_back=0), verbosity=0)

Parameters for an SVD-based DMRG compression algorithm: singular values are truncated according to the truncation scheme `trunc` (any `TruncationScheme`),
while `verbosity` controls the verbosity of the output. The scheme is accessible through the `trunc` field.
"""
struct SVDCompression{T<:TruncationScheme} <: DMRGAlgorithm
	trunc::T
	verbosity::Int
end

"""
    SVDCompression(trunc::TruncationScheme; verbosity=0)

Construct an `SVDCompression` from a `TruncationScheme` (e.g., `truncdimcutoff(D, ϵ)`, `truncdim(D)`, `trunccutoff(ϵ)` or `NoTruncation()`).
"""
SVDCompression(trunc::TruncationScheme; verbosity::Int=0) = SVDCompression(trunc, verbosity)

"""
    SVDCompression(; trunc=truncdimcutoff(D=Defaults.D, ϵ=Defaults.tol, add_back=0), verbosity=0)

Keyword constructor for `SVDCompression`; the default truncation scheme is `truncdimcutoff(D=Defaults.D, ϵ=Defaults.tol, add_back=0)`.
"""
SVDCompression(; trunc::TruncationScheme=truncdimcutoff(D=Defaults.D, ϵ=Defaults.tol, add_back=0), verbosity::Int=0) = SVDCompression(trunc, verbosity)

Base.similar(x::SVDCompression; trunc::TruncationScheme=x.trunc, verbosity::Int=x.verbosity) = SVDCompression(trunc; verbosity=verbosity)

# orthogonalize mps to be left-canonical or right-canonical
abstract type MatrixProductOrthogonalAlgorithm end

"""
	Orthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme}

Configuration of the orthogonalization scheme, used by orthogonalization algorithms such as `leftorth!`, `rightorth!`, and `canonicalize!`.

# Fields
- `orth::A`: underlying orthogonalization algorithm (`QR` or `SVD`)
- `trunc::T`: truncation scheme (`TruncationScheme`); only effective with `SVD`, truncation has no effect with `QR`
- `normalize::Bool`: whether to normalize
- `verbosity::Int`: verbosity level

Main constructor: `Orthogonalize(; alg=SVD(), trunc=NoTruncation(), normalize=false, verbosity=0)`.
"""
struct Orthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme} <: MatrixProductOrthogonalAlgorithm
	orth::A
	trunc::T
	normalize::Bool
	verbosity::Int
end
Orthogonalize(a::Union{QR, SVD}, trunc::TruncationScheme; normalize::Bool=false, verbosity::Int=0) = Orthogonalize(a, trunc, normalize, verbosity)
Orthogonalize(a::Union{QR, SVD}; trunc::TruncationScheme=NoTruncation(), normalize::Bool=false, verbosity::Int=0) = Orthogonalize(a, trunc, normalize, verbosity)
Orthogonalize(; alg::Union{QR, SVD} = SVD(), trunc::TruncationScheme=NoTruncation(), normalize::Bool=false, verbosity::Int=0) = Orthogonalize(alg, trunc, normalize, verbosity)

const AllowedInitGuesses = (:svd, :pre, :rand)

# truncation schemes carrying an explicit maximum bond dimension `D`; `DMRG1`
# requires one of these, since `D` seeds the initial guess of the sweeps
const TruncationWithD = Union{TruncateDim, TruncateDimCutoff}

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
- `trunc::TruncationWithD`: truncation scheme carrying a maximum bond dimension `D` (constructed with `truncdimcutoff(D, ϵ)` or `truncdim(D)`); schemes without a dimension cap (`trunccutoff`, `NoTruncation`) are not allowed
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
