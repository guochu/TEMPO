abstract type AbstractPTLattice{O <: FockOrdering} end

index(x::AbstractPTLattice, args...; kwargs...) = error("index not implemented for pt lattice type $(typeof(x))")
"""
    OrderingStyle(::Type{<:AbstractPTLattice{O}}) where O

Return the Fock ordering instance `O()` associated with the PT lattice.
"""
OrderingStyle(::Type{<:AbstractPTLattice{O}}) where O = O()
"""
    LayoutStyle(::Type{<:AbstractPTLattice{O}}) where O

Return the layout style of the PT lattice, determined by the Fock ordering `O` (e.g. [`TimeLocalLayout`](@ref)).
"""
LayoutStyle(::Type{<:AbstractPTLattice{O}}) where O = LayoutStyle(O)
"""
    ImaginaryTimeOrderingStyle(::Type{<:AbstractPTLattice{O}}) where O

Return the time-ordering style of the PT lattice on the imaginary-time branch, determined by the Fock ordering `O`.
"""
ImaginaryTimeOrderingStyle(::Type{<:AbstractPTLattice{O}}) where O = ImaginaryTimeOrderingStyle(O)
"""
    RealTimeOrderingStyle(::Type{<:AbstractPTLattice{O}}) where O

Return the time-ordering style of the PT lattice on the real-time branch, determined by the Fock ordering `O`.
"""
RealTimeOrderingStyle(::Type{<:AbstractPTLattice{O}}) where O = RealTimeOrderingStyle(O)
OrderingStyle(x::AbstractPTLattice) = OrderingStyle(typeof(x))
LayoutStyle(x::AbstractPTLattice) = LayoutStyle(typeof(x))
ImaginaryTimeOrderingStyle(x::AbstractPTLattice) = ImaginaryTimeOrderingStyle(typeof(x))
RealTimeOrderingStyle(x::AbstractPTLattice) = RealTimeOrderingStyle(typeof(x))
"""
    branches(lattice::AbstractPTLattice)

Return the tuple of branch symbols of the PT lattice: `(:τ,)` for the imaginary-time contour, `(:+, :-)` for the real-time contour, and `(:+, :-, :τ)` for the mixed contour.
"""
branches(lattice::AbstractPTLattice) = branches(typeof(lattice))
phydim(lattice::AbstractPTLattice) = lattice.d

include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


"""
    PTLattice(; contour::Symbol, kwargs...)

Construct a process tensor (PT) lattice according to the contour type.

# Arguments
- `contour::Symbol`: contour type, either `:real` (equivalent to `:Keldysh`), `:imag`, or `:mixed` (equivalent to `:Kadanoff`).
- Remaining keyword arguments are forwarded to the lattice constructor of the corresponding contour.

# Returns
- [`RealPTLattice`](@ref) if `contour == :real`;
- [`ImagPTLattice`](@ref) if `contour == :imag`;
- [`MixedPTLattice`](@ref) if `contour == :mixed`.
"""
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

"""
    vacuumstate(::Type{T}, x::AbstractPTLattice) where {T<:Number}

Construct a vacuum PT state of element type `T` with the same length as the lattice `x`.
"""
vacuumstate(::Type{T}, x::AbstractPTLattice) where {T<:Number} = ProcessTensor(T, length(x), d=x.d)
"""
    vacuumstate(x::AbstractPTLattice)

Construct a vacuum PT state with the same scalar type as the lattice `x`.
"""
vacuumstate(x::AbstractPTLattice) = vacuumstate(scalartype(x), x)

include("integrate.jl")