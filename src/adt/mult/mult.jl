include("svdmult.jl")
include("iterativemult.jl")

mult(x::ADT, y::ADT, alg::SVDCompression) = mult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
mult(x::ADT, y::ADT, alg::DMRGMultAlgorithm) = iterativemult(x, y, alg)


mult!(x::ADT, y::ADT, alg::SVDCompression) = mult!(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
# mult!(x::ADT, y::ADT, alg::DMRGMultAlgorithm) = iterativemult(x, y, alg)

# const DefaultMultAlg = DMRGMult1(DefaultITruncation)
const DefaultMultAlg = SVDCompression(DefaultITruncation)