# abstract type PartialIFStyle end

# struct TriangularPartialIF <: PartialIFStyle end
# struct RectangularPartialIF <: PartialIFStyle end

include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


hybriddynamics(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr, bs; kwargs...)
function hybriddynamics(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(bs), scalartype(corr))
	return hybriddynamics!(vacuumstate(lattice), lattice, corr, bs; kwargs...)
end 


hybriddynamics_naive(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...) = hybriddynamics_naive!(copy(gmps), lattice, corr, bs; kwargs...)
function hybriddynamics_naive(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(bs), scalartype(corr))
	hybriddynamics_naive!(vacuumstate(T, lattice), lattice, corr, bs; kwargs...)
end 



# naive implementation with N^2 gate operations
function hybriddynamics_naive!(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb; 
								trunc::TruncationScheme=DefaultITruncation)
	for b1 in branches(lattice)
		k1 = (b1 == :τ) ? lattice.Nτ : lattice.Nt
		for i in 1:k1
			ind1 = ContourIndex(i, branch=b1)
			tmp = partialif_naive(lattice, ind1, corr, hyb, trunc=trunc)
			gmps = mult!(gmps, tmp, trunc=trunc)			
		end
	end
	return gmps
end
