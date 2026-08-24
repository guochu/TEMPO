include("svdmult.jl")
include("iterativemult.jl")


"""
	mult(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression)

Compute the product of two MPO using the `SVDCompression` compression scheme.

Equivalent to `mult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)`.

# Returns
The product `ProcessTensor`.
"""
mult(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression) = mult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
"""
	mult(x::ProcessTensor, y::ProcessTensor, alg::DMRGMultAlgorithm)

Compute the product of two MPO using a DMRG iterative algorithm (e.g., `DMRGMult1`).

# Returns
The product `ProcessTensor`.
"""
mult(x::ProcessTensor, y::ProcessTensor, alg::DMRGMultAlgorithm) = iterativemult(x, y, alg)


"""
	mult!(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression)

Compute the product of two MPO in place using the `SVDCompression` scheme, storing the result in `x`.

Equivalent to `mult!(x, y, trunc=alg.trunc, verbosity=alg.verbosity)`.

# Returns
`x` itself.
"""
mult!(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression) = mult!(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
"""
	mult!(x::ProcessTensor, y::ProcessTensor, alg::DMRGMultAlgorithm)

Compute the product of two MPO in place using a DMRG iterative algorithm (e.g., `DMRGMult1`), storing the result in `x`.

# Returns
`x` itself.
"""
function mult!(x::ProcessTensor, y::ProcessTensor, alg::DMRGMultAlgorithm)
	r = iterativemult(x, y, alg)
	return copy!(x, r)
end 
