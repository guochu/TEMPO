abstract type MPSAlgorithm end
abstract type DMRGAlgorithm <: MPSAlgorithm end



"""
    SVDCompression(D, tol, verbosity=0)
    SVDCompression(; D=Defaults.D, tol=Defaults.tol, verbosity=0)
    SVDCompression(trunc::TruncationDimCutoff; verbosity=0)

Parameters for an SVD-based DMRG compression algorithm: singular values are truncated by the maximum dimension `D` and the truncation error `tol`,
while `verbosity` controls the verbosity of the output. The equivalent `TruncationDimCutoff` is accessible through the `trunc` property.
"""
struct SVDCompression <: DMRGAlgorithm
	D::Int 
	tol::Float64 
	verbosity::Int 
end

"""
    SVDCompression(; D=Defaults.D, tol=Defaults.tol, verbosity=0)

Keyword constructor for `SVDCompression`; the default dimension and error are taken from `Defaults`.
"""
SVDCompression(; D::Int=Defaults.D, tol::Real=Defaults.tol, verbosity::Int=0) = SVDCompression(D, convert(Float64, tol), verbosity)

"""
    SVDCompression(trunc::TruncationDimCutoff; verbosity=0)

Construct an `SVDCompression` from a `TruncationDimCutoff`; the dimension and error are taken from `trunc`.
"""
SVDCompression(trunc::TruncationDimCutoff; verbosity::Int=0) = SVDCompression(D=trunc.D, tol=trunc.ϵ, verbosity=verbosity)
Base.similar(x::SVDCompression; D::Int=x.D, tol::Float64=x.tol, verbosity::Int=x.verbosity) = SVDCompression(D=D, tol=tol, verbosity=verbosity)

function Base.getproperty(x::SVDCompression, s::Symbol)
	if s == :trunc
		return get_trunc(x)
	elseif s == :ϵ
		return x.tol
	else
		getfield(x, s)
	end
end

get_trunc(alg::SVDCompression) = truncdimcutoff(D=alg.D, ϵ=alg.tol, add_back=0)

# compress!(h::MPO, alg::SVDCompression) = canonicalize!(h, alg=Orthogonalize(SVD(), get_trunc(alg); normalize=false))
# compress!(h::MPO, alg::Deparallelise) = deparallel!(h, tol=alg.tol, verbosity=alg.verbosity)
# compress!(h::MPO; alg::DMRGAlgorithm = Deparallelise()) = compress!(h, alg)
# compress!(psi::MPS, alg::SVDCompression) = canonicalize!(psi, alg=Orthogonalize(trunc=get_trunc(alg), normalize=false))

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

"""
	DMRG1 <: DMRGAlgorithm

Configuration of an MPS/MPO product compression algorithm based on DMRG iterative sweeping.

# Fields
- `trunc::TruncationDimCutoff`: truncation scheme (bond dimension and truncation error)
- `maxiter::Int`: maximum number of iterations
- `tol::Float64`: convergence tolerance
- `initguess::Symbol`: initial guess, one of `:svd`, `:pre`, `:rand`
- `verbosity::Int`: verbosity level
- `callback::Function`: callback function

Main constructor: `DMRG1(trunc; maxiter=5, tol=1e-12, initguess=:svd, verbosity=0, callback=Returns(nothing))`.
"""
struct DMRG1 <: DMRGAlgorithm
	trunc::TruncationDimCutoff
	maxiter::Int
	tol::Float64
	initguess::Symbol
	verbosity::Int
	callback::Function
end
"""
	DMRG1(trunc::TruncationDimCutoff; maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0, callback::Function=Returns(nothing))

Construct a `DMRG1` algorithm configuration.

# Arguments
- `trunc::TruncationDimCutoff`: truncation scheme (can be constructed with `truncdimcutoff(D, ϵ)`)
- `maxiter::Int`: maximum number of iterations
- `tol::Float64`: convergence tolerance
- `initguess::Symbol`: initial guess, must be one of `:svd`, `:pre`, `:rand`, otherwise an `ArgumentError` is thrown
- `verbosity::Int`: verbosity level
- `callback::Function`: callback function
"""
function DMRG1(trunc::TruncationDimCutoff; maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0, callback::Function=Returns(nothing))
	(initguess in AllowedInitGuesses) || throw(ArgumentError("initguess must be one of $(AllowedInitGuesses)"))
	return DMRG1(trunc, maxiter, tol, initguess, verbosity, callback)
end
"""
	DMRG1(; trunc::TruncationDimCutoff=DefaultITruncation, kwargs...)

Construct a `DMRG1` from keyword arguments, with default truncation scheme `DefaultITruncation`.
"""
DMRG1(; trunc::TruncationDimCutoff=DefaultITruncation, kwargs...) = DMRG1(trunc; kwargs...)
Base.similar(x::DMRG1; trunc::TruncationDimCutoff=x.trunc, maxiter::Int=x.maxiter, tol::Float64=x.tol, initguess::Symbol=x.initguess, verbosity::Int=x.verbosity, callback=x.callback) = DMRG1(
			trunc=trunc, maxiter=maxiter, tol=tol, initguess=initguess, verbosity=verbosity, callback=callback)

function Base.getproperty(x::DMRGAlgorithm, s::Symbol)
	if s == :D
		return x.trunc.D
	elseif s == :ϵ
		return x.trunc.ϵ
	else
		getfield(x, s)
	end
end
