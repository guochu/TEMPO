"""
    ImagADTLattice{O<:ImagFockOrdering} <: AbstractADTLattice{O}

Abstract ADT lattice type on the imaginary-time (Matsubara) contour: scalar type `Float64`, containing only the `:τ` branch.
"""
abstract type ImagADTLattice{O<:ImagFockOrdering} <: AbstractADTLattice{O} end
TO.scalartype(::Type{<:ImagADTLattice}) = Float64
"""
    branches(::Type{<:ImagADTLattice})

Tuple of branch symbols of the imaginary-time ADT lattice, `(:τ,)`.
"""
branches(::Type{<:ImagADTLattice}) = (:τ,)
TimeOrderingStyle(x::ImagADTLattice) = ImaginaryTimeOrderingStyle(x)

"""
    ImagADTLattice1Order{O<:ImagFockOrdering} <: ImagADTLattice{O}

First-order ADT lattice on the imaginary-time contour: discretizes the imaginary-time interval `[0, β]` into `N` steps of size `δτ = β/N`.

# Fields
- `δτ::Float64`: imaginary-time step size.
- `d::Int`: physical dimension.
- `N::Int`: number of imaginary-time discretization steps.
- `ordering::O`: Fock ordering (default [`M2M1`](@ref)).

# Derived attributes
- `β = N * δτ`: inverse temperature; `T = 1/β`.
- `τs = 0:δτ:β`: imaginary-time grid.
- `k = kτ = N + 1`: number of lattice indices (including boundary points).
"""
struct ImagADTLattice1Order{O<:ImagFockOrdering} <: ImagADTLattice{O}
	δτ::Float64
	d::Int
	N::Int
	ordering::O

	ImagADTLattice1Order(δτ::Real, d::Int, N::Int, ordering::ImagFockOrdering) = new{typeof(ordering)}(float(δτ), d, N, ordering)
end

"""
    ImagADTLattice1Order(; δτ::Real, N::Int, d::Int=2, ordering::ImagFockOrdering=M2M1())

Convenience constructor for a first-order imaginary-time ADT lattice.
"""
ImagADTLattice1Order(; δτ::Real, N::Int, d::Int=2, ordering::ImagFockOrdering=M2M1()) = ImagADTLattice1Order(δτ, d, N, ordering)
Base.similar(x::ImagADTLattice1Order; δτ::Real=x.δτ, d::Int=x.d, N::Int=x.N, ordering::ImagFockOrdering=x.ordering) = ImagADTLattice1Order(δτ, d, N, ordering)

"""
    ImagADTLattice(; order::Int=1, kwargs...)

Construct an imaginary-time ADT lattice; `order` specifies the splitting order (currently only first order is supported, i.e. `order == 1`).
"""
function ImagADTLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return ImagADTLattice1Order(; kwargs...)
	else
		error("Second orderr ImagGrassmannLattice not implemented")
	end
end

Base.length(x::ImagADTLattice) = x.k
function Base.getproperty(x::ImagADTLattice1Order, s::Symbol)
	if s == :τs
		return 0:x.δτ:x.N*x.δτ
	elseif s == :β
		return x.N * x.δτ
	elseif s == :T
		return 1 / x.β
	elseif s == :Nτ
		return x.N
	elseif (s == :k) || ( s == :kτ)
		return x.N + 1
	else
		getfield(x, s)
	end
end


function index(x::ImagADTLattice{<:M2M1}, i::Int; branch::Symbol=:τ)
	@boundscheck begin
		(1 <= i <= x.k) || throw(BoundsError(1:x.k, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	return x.k-i+1
end



"""
    indexmappings(lattice::ImagADTLattice)

Return a `Dict{Tuple{Int, Symbol}, Int}` mapping `(time step, branch)` to the global index of the ADT lattice.
"""
function indexmappings(lattice::ImagADTLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.k
		r[(i, :τ)] = index(lattice, i)
	end
	return r
end
