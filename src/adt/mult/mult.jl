include("svdmult.jl")
include("iterativemult.jl")

"""
	mult(x::ADT, y::ADT, alg::SVDCompression)

Compute the product of two MPS using the `SVDCompression` compression scheme.

Equivalent to `mult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)`.

# Returns
The product `ADT`.
"""
mult(x::ADT, y::ADT, alg::SVDCompression) = mult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
"""
	mult(x::ADT, y::ADT, alg::DMRGMultAlgorithm)

Compute the product of two MPS using a DMRG iterative algorithm (e.g., `DMRGMult1`).

# Returns
The product `ADT`.
"""
mult(x::ADT, y::ADT, alg::DMRGMultAlgorithm) = iterativemult(x, y, alg)


"""
	mult!(x::ADT, y::ADT, alg::SVDCompression)

Compute the product of two MPS in place using the `SVDCompression` scheme, storing the result in `x`.

Equivalent to `mult!(x, y, trunc=alg.trunc, verbosity=alg.verbosity)`.

# Returns
`x` itself.
"""
mult!(x::ADT, y::ADT, alg::SVDCompression) = mult!(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
"""
	mult!(x::ADT, y::ADT, alg::DMRGMultAlgorithm)

Compute the product of two MPS in place using a DMRG iterative algorithm (e.g., `DMRGMult1`), storing the result in `x`.

# Returns
`x` itself.
"""
function mult!(x::ADT, y::ADT, alg::DMRGMultAlgorithm)
	r = iterativemult(x, y, alg)
	return copy!(x, r)
end 

const DefaultMultAlg = DMRGMult1(DefaultITruncation)
# const DefaultMultAlg = SVDCompression(DefaultITruncation)