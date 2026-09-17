
#default settings
module Defaults
	const maxiter = 100 # for DMRG iteration
	const D = 100 # default bond dimension 
	const tolgauge = 1e-14 # for MPS truncation
	const tol = 1e-12 # for DMRG iteration
	const tollanczos = 1.0e-10 # for lanczos eigensolver
	const tolexp = 1.0e-8 # for local eigen in DMRG
	const verbosity = 1
end

const DefaultITruncation = truncdimcutoff(D=Defaults.D, ϵ=Defaults.tolgauge, add_back=0) # for IF construction and general MPS/MPO compression
const DefaultKTruncation = trunccutoff(Defaults.tolgauge) # system dynamics, initial-state absorption and MPO compression

const DefaultMultAlg = DMRG1(DefaultITruncation) # default compression algorithm for MPS/MPO multiplication
# const DefaultMultAlg = SVDCompression(DefaultITruncation)
