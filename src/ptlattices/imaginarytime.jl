"""
    AbstractImagPTLattice{O<:ImagFockOrdering} <: AbstractPTLattice{O}

Abstract PT lattice type on the imaginary-time (Matsubara) contour: scalar type `Float64`, containing only the `:τ` branch.
"""
abstract type AbstractImagPTLattice{O<:ImagFockOrdering} <: AbstractPTLattice{O} end
TO.scalartype(::Type{<:AbstractImagPTLattice}) = Float64
"""
    branches(::Type{<:AbstractImagPTLattice})

Tuple of branch symbols of the imaginary-time PT lattice, `(:τ,)`.
"""
branches(::Type{<:AbstractImagPTLattice}) = (:τ,)
TimeOrderingStyle(x::AbstractImagPTLattice) = ImaginaryTimeOrderingStyle(x)

"""
    ImagPTLattice{O<:ImagFockOrdering} <: AbstractImagPTLattice{O}

PT lattice on the imaginary-time (Matsubara) contour: scalar type `Float64`, containing only the `:τ` branch.
The lattice discretizes the imaginary-time interval `[0, β]` into `N` steps of size `δτ = β/N`.

# Fields
- `δτ::Float64`: imaginary-time step size.
- `d::Int`: physical dimension.
- `N::Int`: number of imaginary-time discretization steps.
- `ordering::O`: Fock ordering (default [`M2M1`](@ref)).

# Derived attributes
- `β = N * δτ`: inverse temperature; `T = 1/β`.
- `τs = 0:δτ:β`: imaginary-time grid.
- `Nτ = N`: number of imaginary-time discretization steps (including boundary points `N + 1`).
"""
struct ImagPTLattice{O<:ImagFockOrdering} <: AbstractImagPTLattice{O}
	δτ::Float64
	d::Int
	N::Int
	ordering::O

	ImagPTLattice(δτ::Real, d::Int, N::Int, ordering::ImagFockOrdering) = new{typeof(ordering)}(float(δτ), d, N, ordering)
end

"""
    ImagPTLattice(; δτ::Real, N::Int, d::Int=2, ordering::ImagFockOrdering=M2M1())

Convenience constructor for an imaginary-time PT lattice.
"""
ImagPTLattice(; δτ::Real, N::Int, d::Int=2, ordering::ImagFockOrdering=M2M1()) = ImagPTLattice(δτ, d, N, ordering)
Base.similar(x::ImagPTLattice; δτ::Real=x.δτ, d::Int=x.d, N::Int=x.N, ordering::ImagFockOrdering=x.ordering) = ImagPTLattice(δτ, d, N, ordering)

Base.length(x::ImagPTLattice) = x.N
function Base.getproperty(x::ImagPTLattice, s::Symbol)
	if s == :τs
		return 0:x.δτ:x.N*x.δτ
	elseif s == :β
		return x.N * x.δτ
	elseif s == :T
		return 1 / x.β
	elseif s == :Nτ
		return x.N
	else
		getfield(x, s)
	end
end


"""
    getindex(x::ImagPTLattice{<:M2M1}, a::BosonicContourIndex)

Map the contour index `a` (time step `a.j`, branch `:τ`) to the global lattice position.
"""
function Base.getindex(x::ImagPTLattice{<:M2M1}, a::BosonicContourIndex)
	@boundscheck begin
		(1 <= a.j <= x.N) || throw(BoundsError(1:x.N, a.j))
		(branch(a) == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	return x.N-a.j+1
end



"""
    indexmappings(lattice::ImagPTLattice)

Return a `Dict{Tuple{Int, Symbol}, Int}` mapping `(time step, branch)` to the global index of the PT lattice.
"""
function indexmappings(lattice::ImagPTLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.N
		r[(i, :τ)] = lattice[BosonicContourIndex(i, :τ)]
	end
	return r
end
