abstract type AbstractPTLattice{O <: FockOrdering} end

index(x::AbstractPTLattice, args...; kwargs...) = error("index not implemented for pt lattice type $(typeof(x))")
OrderingStyle(::Type{<:AbstractPTLattice{O}}) where O = O()
LayoutStyle(::Type{<:AbstractPTLattice{O}}) where O = LayoutStyle(O)
ImaginaryTimeOrderingStyle(::Type{<:AbstractPTLattice{O}}) where O = ImaginaryTimeOrderingStyle(O)
RealTimeOrderingStyle(::Type{<:AbstractPTLattice{O}}) where O = RealTimeOrderingStyle(O)
OrderingStyle(x::AbstractPTLattice) = OrderingStyle(typeof(x))
LayoutStyle(x::AbstractPTLattice) = LayoutStyle(typeof(x))
ImaginaryTimeOrderingStyle(x::AbstractPTLattice) = ImaginaryTimeOrderingStyle(typeof(x))
RealTimeOrderingStyle(x::AbstractPTLattice) = RealTimeOrderingStyle(typeof(x))
branches(lattice::AbstractPTLattice) = branches(typeof(lattice))
phydim(lattice::AbstractPTLattice) = lattice.d

include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


function PTLattice(; contour::Symbol, kwargs...)
	(contour in (:real, :imag, :Keldysh, :mixed, :Kadanoff)) || throw(ArgumentError("contour must be :real (equivalentlt :Keldysh), :imag or :mixed (equivalentlt :KadanoffBaym)"))
	if (contour == :real) || (contour == :Keldysh)
		return RealPTLattice(; kwargs...)
	elseif contour == :imag
		return ImagPTLattice(; kwargs...)
	else
		return MixedPTLattice(; kwargs...)
	end
end

Base.getindex(lat::AbstractPTLattice, a::ContourIndex) = index(lat, a.j, branch=branch(a))

vacuumstate(x::AbstractPTLattice) = ProcessTensor(scalartype(x), length(x), d=x.d)


include("integrate.jl")