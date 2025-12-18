# abstract type PartialIFStyle end

# struct TriangularPartialIF <: PartialIFStyle end
# struct RectangularPartialIF <: PartialIFStyle end


include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


hybriddynamics(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::HybridizationStyle; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr, bs; kwargs...)
hybriddynamics(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::HybridizationStyle; kwargs...) = hybriddynamics!(vacuumstate(lattice), lattice, corr, bs; kwargs...)


hybriddynamics_naive(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::HybridizationStyle; kwargs...) = hybriddynamics_naive!(copy(gmps), lattice, corr, bs; kwargs...)
hybriddynamics_naive(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::HybridizationStyle; kwargs...) = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr, bs; kwargs...)
