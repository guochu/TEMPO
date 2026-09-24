include("iterativemult.jl")

"""
	mult(x::ADT, y::ADT; trunc::TruncationScheme=DefaultITruncation, verbosity::Int=0)

Compute the (compressed) product of two MPS and return a new `ADT` (inputs are not modified).

The SVD compression route delegates to FiniteMPSAlgorithms' `hadamard`: the exact
(Hadamard) product `x ⊙ y` is formed and compressed by a single SVD sweep under `trunc`.

# Returns
The product `ADT`.
"""
mult(x::ADT, y::ADT; trunc::TruncationScheme = DefaultITruncation, verbosity::Int = 0) =
	mult(x, y, SVDCompression(trunc; verbosity))
"""
	mult(x::ADT, y::ADT, alg::SVDCompression)

Compute the product of two MPS using the `SVDCompression` compression scheme;
delegates to FiniteMPSAlgorithms' `hadamard` on the inner `CanonicalMPS` payloads.
"""
mult(x::ADT, y::ADT, alg::SVDCompression) = ADT(hadamard(x.parent, y.parent, alg))
"""
	mult(x::ADT, y::ADT, alg::DMRGAlgorithm)

Compute the product of two MPS using a DMRG iterative algorithm (e.g., `DMRG1`).
"""
mult(x::ADT, y::ADT, alg::DMRGAlgorithm) = iterativemult(x, y, alg)

"""
	mult!(x::ADT, y::ADT; trunc::TruncationScheme=DefaultITruncation, verbosity::Int=0)

Compute the (compressed) product of two MPS in place on `x`, and return `x`.
"""
mult!(x::ADT, y::ADT; trunc::TruncationScheme = DefaultITruncation, verbosity::Int = 0) =
	copy!(x, mult(x, y; trunc, verbosity))
"""
	mult!(x::ADT, y::ADT, alg::SVDCompression)

Compute the product of two MPS in place using the `SVDCompression` scheme, storing the result in `x`.
"""
mult!(x::ADT, y::ADT, alg::SVDCompression) = copy!(x, mult(x, y, alg))
"""
	mult!(x::ADT, y::ADT, alg::DMRGAlgorithm)

Compute the product of two MPS in place using a DMRG iterative algorithm (e.g., `DMRG1`), storing the result in `x`.
"""
function mult!(x::ADT, y::ADT, alg::DMRGAlgorithm)
	r = iterativemult(x, y, alg)
	return copy!(x, r)
end
