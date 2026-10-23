"""
    AbstractRealPTLattice{O<:RealFockOrdering} <: AbstractPTLattice{O}

Abstract PT lattice type on the real-time (Keldysh) contour: scalar type `ComplexF64`, containing the `:+` and `:-` branches.
"""
abstract type AbstractRealPTLattice{O<:RealFockOrdering} <: AbstractPTLattice{O} end
TO.scalartype(::Type{<:AbstractRealPTLattice}) = ComplexF64
"""
    branches(::Type{<:AbstractRealPTLattice})

Tuple of branch symbols of the real-time PT lattice, `(:+, :-)`.
"""
branches(::Type{<:AbstractRealPTLattice}) = (:+, :-)
TimeOrderingStyle(x::AbstractRealPTLattice) = RealTimeOrderingStyle(x)

"""
    RealPTLattice{O<:RealFockOrdering} <: AbstractRealPTLattice{O}

PT lattice on the real-time (Keldysh) contour: scalar type `ComplexF64`, containing the `:+` and `:-` branches.
The lattice discretizes the real-time interval `[0, t]` into `N` steps of size `δt = t/N`.

# Fields
- `δt::Float64`: real-time step size.
- `d::Int`: physical dimension.
- `N::Int`: number of real-time discretization steps.
- `ordering::O`: Fock ordering (default [`M2m2M1m1`](@ref)).

# Derived attributes
- `t = N * δt`: total evolution time; `ts = 0:δt:t`: real-time grid.
- `Nt = N`: number of real-time discretization steps.
"""
struct RealPTLattice{O<:RealFockOrdering} <: AbstractRealPTLattice{O}
	δt::Float64
	d::Int
	N::Int
	ordering::O

	RealPTLattice(δt::Real, d::Int, N::Int, ordering::RealFockOrdering) = new{typeof(ordering)}(float(δt), d, N, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
"""
    RealPTLattice(; δt::Real, N::Int, d::Int=2, ordering::RealFockOrdering=M2m2M1m1())

Convenience constructor for a real-time PT lattice.
"""
RealPTLattice(; δt::Real, N::Int, d::Int=2, ordering::RealFockOrdering=M2m2M1m1()) = RealPTLattice(δt, d, N, ordering)
Base.similar(x::RealPTLattice; δt::Real=x.δt, d::Int=x.d, N::Int=x.N, ordering::RealFockOrdering=x.ordering) = RealPTLattice(δt, d, N, ordering)

# k is the number of discretization, nbands is the number of bands
# pos is the position within a band
# k+1 due to the grassmann number on the boundary for the final trace

function Base.getproperty(x::RealPTLattice, s::Symbol)
	if s == :t
		return x.N * x.δt
	elseif s == :ts
		return 0:x.δt:x.N*x.δt
	elseif s == :Nt
		return x.N
	else
		getfield(x, s)
	end
end
Base.length(x::RealPTLattice) = 2*x.N


"""
    getindex(x::RealPTLattice{<:M2m2M1m1}, a::BosonicContourIndex)

Map the contour index `a` (time step `a.j`, branch `:+` or `:-`) to the global lattice position.
"""
function Base.getindex(x::RealPTLattice{<:M2m2M1m1}, a::BosonicContourIndex)
	@boundscheck begin
		(1 <= a.j <= x.N) || throw(BoundsError(1:x.N, a.j))
		(branch(a) in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	ifelse(branch(a) == :+, TL-2a.j+1, TL-2a.j+2)
end


# key is timestep, conj, branch, band
"""
    indexmappings(lattice::RealPTLattice)

Return a `Dict{Tuple{Int, Symbol}, Int}` mapping `(time step, branch)` to the global index of the PT lattice.
"""
function indexmappings(lattice::RealPTLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.N
		for f in (:+, :-)
			r[(i, f)] = lattice[BosonicContourIndex(i, f)]
		end
	end
	return r
end

