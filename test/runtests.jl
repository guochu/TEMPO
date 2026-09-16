# Test suite organized in three parts (following the GTEMPO layout):
#   (i)   api/       : API unit tests of all exported functions (no physical references)
#   (ii)  adtmodels/ : functional tests of ADT-lattice models against ED / analytic references
#   (iii) ptmodels/  : functional tests of PT-lattice models against ED / analytic references

using Test, LinearAlgebra, TensorOperations, ImpurityModelBase, QuAPI
using Random

using TEMPO
using TEMPO: QR, QRpos, LQ, LQpos, SVD, SDD, Polar, TimeAscending, TimeDscending, Zvalue2

Random.seed!(20260916)

include("util.jl")

@testset "TEMPO" begin
	@testset "api" begin
		include("api/api.jl")
	end
	@testset "adtmodels" begin
		include("adtmodels/adtmodels.jl")
	end
	@testset "ptmodels" begin
		include("ptmodels/ptmodels.jl")
	end
end
