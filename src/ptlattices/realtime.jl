"""
    RealPTLattice{O<:RealFockOrdering} <: AbstractPTLattice{O}

Abstract PT lattice type on the real-time (Keldysh) contour: scalar type `ComplexF64`, containing the `:+` and `:-` branches.
"""
abstract type RealPTLattice{O<:RealFockOrdering} <: AbstractPTLattice{O} end
TO.scalartype(::Type{<:RealPTLattice}) = ComplexF64
"""
    branches(::Type{<:RealPTLattice})

Tuple of branch symbols of the real-time PT lattice, `(:+, :-)`.
"""
branches(::Type{<:RealPTLattice}) = (:+, :-)
TimeOrderingStyle(x::RealPTLattice) = RealTimeOrderingStyle(x)

# k is the number of discretization, nbands is the number of bands
# pos is the position within a band
# k+1 due to the grassmann number on the boundary for the final trace

"""
	struct RealPTLattice1Order{O<:RealFockOrdering} <: RealPTLattice{O}

First-order PT lattice on the real-time contour: discretizes the real-time interval `[0, t]` into `N` steps of size `δt = t/N`.

# Fields
- `δt::Float64`: real-time step size.
- `d::Int`: physical dimension.
- `N::Int`: number of real-time discretization steps.
- `ordering::O`: Fock ordering (default [`M2m2M1m1`](@ref)).

# Derived attributes
- `t = N * δt`: total evolution time; `ts = 0:δt:t`: real-time grid.
- `Nt = N`: number of real-time discretization steps.
"""
struct RealPTLattice1Order{O<:RealFockOrdering} <: RealPTLattice{O}
	δt::Float64
	d::Int
	N::Int
	ordering::O

	RealPTLattice1Order(δt::Real, d::Int, N::Int, ordering::RealFockOrdering) = new{typeof(ordering)}(float(δt), d, N, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
"""
    RealPTLattice1Order(; δt::Real, N::Int, d::Int=2, ordering::RealFockOrdering=M2m2M1m1())

Convenience constructor for a first-order real-time PT lattice.
"""
RealPTLattice1Order(; δt::Real, N::Int, d::Int=2, ordering::RealFockOrdering=M2m2M1m1()) = RealPTLattice1Order(δt, d, N, ordering)
Base.similar(x::RealPTLattice1Order; δt::Real=x.δt, d::Int=x.d, N::Int=x.N, ordering::RealFockOrdering=x.ordering) = RealPTLattice1Order(δt, d, N, ordering)

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
    RealPTLattice(; order::Int=1, kwargs...)

Construct a real-time PT lattice; `order` specifies the splitting order (currently only first order is supported, i.e. `order == 1`).
"""
function RealPTLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return RealPTLattice1Order(; kwargs...)
	else
		error("Second orderr RealPTLattice not implemented")
	end
end


function index(x::RealPTLattice{<:M2m2M1m1}, i::Int; branch::Symbol=:+)
	@boundscheck begin
		(1 <= i <= x.N) || throw(BoundsError(1:x.N, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	ifelse(branch == :+, TL-2i+1, TL-2i+2)
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
			r[(i, f)] = index(lattice, i, branch=f)
		end
	end
	return r
end

