# abstract type PartialIFStyle end

# struct TriangularPartialIF <: PartialIFStyle end
# struct RectangularPartialIF <: PartialIFStyle end

include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


"""
	hybriddynamics(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb; kwargs...)
	hybriddynamics(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb; kwargs...)

Construct the influence functional for the given lattice, bath correlation function and additive (diagonal) system-bath coupling. If no initial `gmps` is supplied, `vacuumstate(lattice)` is used as the initial state.

# Arguments
- `gmps::ADT`: initial augmented density tensor (multiplication is performed on a copy; the original object is left unchanged).
- `lattice::AbstractADTLattice`: ADT contour lattice (imaginary-time/real-time/mixed-time).
- `corr::AbstractCorrelationFunction`: bath correlation function.
- `hyb::AdditiveHyb`: additive (diagonal) system-bath coupling.
- `kwargs...`: keyword arguments passed to `hybriddynamics!` (e.g. `trunc=...`).

# Returns
The influence functional, represented as an `ADT`.

See also [`hybriddynamics!`](@ref) and [`hybriddynamics_naive`](@ref).
"""
hybriddynamics(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr, bs; kwargs...)
function hybriddynamics(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(bs), scalartype(corr))
	return hybriddynamics!(vacuumstate(lattice), lattice, corr, bs; kwargs...)
end 


"""
	hybriddynamics_naive(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb; kwargs...)
	hybriddynamics_naive(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb; kwargs...)

Naive implementation of `hybriddynamics`: construct a partial IF for each contour lattice point and multiply them in successively (N² gate operations in total); simpler but slower than the optimized version.

See also [`hybriddynamics`](@ref) and [`hybriddynamics_naive!`](@ref).
"""
hybriddynamics_naive(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...) = hybriddynamics_naive!(copy(gmps), lattice, corr, bs; kwargs...)
function hybriddynamics_naive(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(bs), scalartype(corr))
	hybriddynamics_naive!(vacuumstate(T, lattice), lattice, corr, bs; kwargs...)
end 



# naive implementation with N^2 gate operations
"""
	hybriddynamics_naive!(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb; trunc::TruncationScheme=DefaultITruncation)

In-place version of `hybriddynamics_naive`: construct a partial IF for each contour lattice point and multiply it directly into `gmps` (N² gate operations in total), truncating after each multiplication according to `trunc`.

# Returns
The modified `gmps`.

See also [`hybriddynamics_naive`](@ref) and [`hybriddynamics`](@ref).
"""
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
