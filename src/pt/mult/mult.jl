include("svdmult.jl")
include("iterativemult.jl")


mult(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression) = mult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
mult(x::ProcessTensor, y::ProcessTensor, alg::DMRGMultAlgorithm) = iterativemult(x, y, alg)


mult!(x::ProcessTensor, y::ProcessTensor, alg::SVDCompression) = mult!(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
function mult!(x::ProcessTensor, y::ProcessTensor, alg::DMRGMultAlgorithm)
	r = iterativemult(x, y, alg)
	return copy!(x, r)
end 
