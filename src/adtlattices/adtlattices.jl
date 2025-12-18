include("fockordering.jl")

abstract type AbstractADTLattice{O <: FockOrdering} end

index(x::AbstractADTLattice, args...; kwargs...) = error("index not implemented for grassmann lattice type $(typeof(x))")
OrderingStyle(::Type{<:AbstractADTLattice{O}}) where O = O()
LayoutStyle(::Type{<:AbstractADTLattice{O}}) where O = LayoutStyle(O)
ImaginaryTimeOrderingStyle(::Type{<:AbstractADTLattice{O}}) where O = ImaginaryTimeOrderingStyle(O)
RealTimeOrderingStyle(::Type{<:AbstractADTLattice{O}}) where O = RealTimeOrderingStyle(O)
OrderingStyle(x::AbstractADTLattice) = OrderingStyle(typeof(x))
LayoutStyle(x::AbstractADTLattice) = LayoutStyle(typeof(x))
ImaginaryTimeOrderingStyle(x::AbstractADTLattice) = ImaginaryTimeOrderingStyle(typeof(x))
RealTimeOrderingStyle(x::AbstractADTLattice) = RealTimeOrderingStyle(typeof(x))
branches(lattice::AbstractADTLattice) = branches(typeof(lattice))
phydim(lattice::AbstractADTLattice) = lattice.d

include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


function ADTLattice(; contour::Symbol, kwargs...)
	(contour in (:real, :imag, :Keldysh, :mixed, :Kadanoff)) || throw(ArgumentError("contour must be :real (equivalentlt :Keldysh), :imag or :mixed (equivalentlt :KadanoffBaym)"))
	if (contour == :real) || (contour == :Keldysh)
		return RealADTLattice(; kwargs...)
	elseif contour == :imag
		return ImagADTLattice(; kwargs...)
	else
		return MixedADTLattice(; kwargs...)
	end
end

Base.getindex(lat::AbstractADTLattice, a::ContourIndex) = index(lat, a.j, branch=branch(a))

vacuumstate(x::AbstractADTLattice) = ADT(scalartype(x), length(x), d=x.d)