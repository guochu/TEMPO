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
