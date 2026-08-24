include("fockordering.jl")

abstract type AbstractADTLattice{O <: FockOrdering} end

index(x::AbstractADTLattice, args...; kwargs...) = error("index not implemented for fock lattice type $(typeof(x))")
"""
    OrderingStyle(::Type{<:AbstractADTLattice{O}}) where O

Return the Fock ordering instance `O()` associated with the ADT lattice.
"""
OrderingStyle(::Type{<:AbstractADTLattice{O}}) where O = O()
"""
    LayoutStyle(::Type{<:AbstractADTLattice{O}}) where O

Return the layout style of the ADT lattice, determined by the Fock ordering `O` (e.g. [`TimeLocalLayout`](@ref)).
"""
LayoutStyle(::Type{<:AbstractADTLattice{O}}) where O = LayoutStyle(O)
"""
    ImaginaryTimeOrderingStyle(::Type{<:AbstractADTLattice{O}}) where O

Return the time-ordering style of the ADT lattice on the imaginary-time branch, determined by the Fock ordering `O`.
"""
ImaginaryTimeOrderingStyle(::Type{<:AbstractADTLattice{O}}) where O = ImaginaryTimeOrderingStyle(O)
"""
    RealTimeOrderingStyle(::Type{<:AbstractADTLattice{O}}) where O

Return the time-ordering style of the ADT lattice on the real-time branch, determined by the Fock ordering `O`.
"""
RealTimeOrderingStyle(::Type{<:AbstractADTLattice{O}}) where O = RealTimeOrderingStyle(O)
OrderingStyle(x::AbstractADTLattice) = OrderingStyle(typeof(x))
LayoutStyle(x::AbstractADTLattice) = LayoutStyle(typeof(x))
ImaginaryTimeOrderingStyle(x::AbstractADTLattice) = ImaginaryTimeOrderingStyle(typeof(x))
RealTimeOrderingStyle(x::AbstractADTLattice) = RealTimeOrderingStyle(typeof(x))
"""
    branches(lattice::AbstractADTLattice)

Return the tuple of branch symbols of the ADT lattice: `(:τ,)` for the imaginary-time contour, `(:+, :-)` for the real-time contour, and `(:+, :-, :τ)` for the mixed contour.
"""
branches(lattice::AbstractADTLattice) = branches(typeof(lattice))
phydim(lattice::AbstractADTLattice) = lattice.d

include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


"""
    ADTLattice(; contour::Symbol, kwargs...)

Construct an augmented density tensor (ADT) lattice according to the contour type.

# Arguments
- `contour::Symbol`: contour type, either `:real` (equivalent to `:Keldysh`), `:imag`, or `:mixed` (equivalent to `:Kadanoff`).
- Remaining keyword arguments are forwarded to the lattice constructor of the corresponding contour (e.g. `δt`, `N`, `δτ`, `d`, etc.).

# Returns
- [`RealADTLattice`](@ref) if `contour == :real`;
- [`ImagADTLattice`](@ref) if `contour == :imag`;
- [`MixedADTLattice`](@ref) if `contour == :mixed`.

# Examples
```julia
lattice = ADTLattice(contour=:real, δt=0.1, N=100)
```
"""
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

"""
    vacuumstate(::Type{T}, x::AbstractADTLattice) where {T<:Number}

Construct a vacuum ADT state of element type `T` with the same length as the lattice `x`.
"""
vacuumstate(::Type{T}, x::AbstractADTLattice) where {T<:Number} = ADT(T, length(x), d=x.d)
"""
    vacuumstate(x::AbstractADTLattice)

Construct a vacuum ADT state with the same scalar type as the lattice `x`.
"""
vacuumstate(x::AbstractADTLattice) = vacuumstate(scalartype(x), x)