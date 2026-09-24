include("iterativemult.jl")

"""
	mult(x::ProcessTensor, y::ProcessTensor; trunc::TruncationScheme=DefaultITruncation, verbosity::Int=0)

Compute the (compressed) product of two MPO and return a new `ProcessTensor` (inputs are not modified).

The SVD compression route delegates to FiniteMPSAlgorithms' `mult`: the exact product
`x.parent * y.parent` is formed and compressed by a single SVD sweep under `trunc`.
"""
mult(x::ProcessTensor, y::ProcessTensor; trunc::TruncationScheme = DefaultITruncation, verbosity::Int = 0) =
	mult(x, y, SVDCompression(trunc; verbosity))
"""
	mult(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression)

Compute the product of two MPO using the `SVDCompression` compression scheme;
delegates to FiniteMPSAlgorithms' `mult` on the inner `CanonicalMPO` payloads.
"""
mult(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression) = ProcessTensor(mult(x.parent, y.parent, alg))
"""
	mult(x::ProcessTensor, y::ProcessTensor, alg::DMRGAlgorithm)

Compute the product of two MPO using a DMRG iterative algorithm (e.g., `DMRG1`).
"""
mult(x::ProcessTensor, y::ProcessTensor, alg::DMRGAlgorithm) = iterativemult(x, y, alg)

"""
	mult!(x::ProcessTensor, y::ProcessTensor; trunc::TruncationScheme=DefaultITruncation, verbosity::Int=0)

Compute the (compressed) product of two MPO in place on `x`, and return `x`.
"""
mult!(x::ProcessTensor, y::ProcessTensor; trunc::TruncationScheme = DefaultITruncation, verbosity::Int = 0) =
	copy!(x, mult(x, y; trunc, verbosity))
"""
	mult!(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression)

Compute the product of two MPO in place using the `SVDCompression` scheme, storing the result in `x`.
"""
mult!(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression) = copy!(x, mult(x, y, alg))
"""
	mult!(x::ProcessTensor, y::ProcessTensor, alg::DMRGAlgorithm)

Compute the product of two MPO in place using a DMRG iterative algorithm (e.g., `DMRG1`), storing the result in `x`.
"""
function mult!(x::ProcessTensor, y::ProcessTensor, alg::DMRGAlgorithm)
	r = iterativemult(x, y, alg)
	return copy!(x, r)
end
