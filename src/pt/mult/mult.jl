include("svdmult.jl")
include("iterativemult.jl")


mult(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression) = mult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)

mult!(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression) = mult!(x, y, trunc=alg.trunc, verbosity=alg.verbosity)