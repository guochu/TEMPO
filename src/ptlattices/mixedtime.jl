"""
    AbstractMixedPTLattice{O<:MixedFockOrdering} <: AbstractPTLattice{O}

Abstract PT lattice type on the mixed (Kadanoff-Baym, real-time + imaginary-time) contour: scalar type `ComplexF64`,
containing the `:+`, `:-` and `:τ` branches.
"""
abstract type AbstractMixedPTLattice{O<:MixedFockOrdering} <: AbstractPTLattice{O} end
TO.scalartype(::Type{<:AbstractMixedPTLattice}) = ComplexF64
"""
    branches(::Type{<:AbstractMixedPTLattice})

Tuple of branch symbols of the mixed PT lattice, `(:+, :-, :τ)`.
"""
branches(::Type{<:AbstractMixedPTLattice}) = (:+, :-, :τ)

"""
    MixedPTLattice{O<:MixedFockOrdering} <: AbstractMixedPTLattice{O}

PT lattice on the mixed (Kadanoff-Baym, real-time + imaginary-time) contour: scalar type `ComplexF64`,
containing the `:+`, `:-` and `:τ` branches. The lattice discretizes the real-time interval `[0, t]` into `Nt` steps
(size `δt`) and the imaginary-time interval `[0, β]` into `Nτ` steps (size `δτ`).

# Fields
- `δt::Float64`: real-time step size.
- `Nt::Int`: number of real-time discretization steps.
- `δτ::Float64`: imaginary-time step size.
- `Nτ::Int`: number of imaginary-time discretization steps.
- `d::Int`: physical dimension.
- `ordering::O`: Fock ordering (default [`M2M1_m1M1m2M2`](@ref)).

# Derived attributes
- `t = Nt * δt`: total evolution time; `β = Nτ * δτ`: inverse temperature; `T = 1/β`.
- `ts = 0:δt:t`: real-time grid; `τs = 0:δτ:β`: imaginary-time grid.

# Indexing convention

Unlike the mixed ADT lattice, the mixed PT lattice stores exactly `Nτ` imaginary-time and `Nt` real-time points per
branch (no extra boundary points). Correlation-function and operator-insertion indices therefore map **directly** to
the lattice indices (`lattice[BosonicContourIndex(i, branch)]` with no shift): the `BosonicContourIndex(i, branch)` of a gate or
inserted operator addresses the lattice point at contour time `(i-1)δ` on the corresponding branch. This is the same
direct mapping as the pure imaginary-time/real-time PT lattices and is used uniformly in `partialif_naive` and the
`ContourOperator` `apply!` on PT lattices.
"""
struct MixedPTLattice{O<:MixedFockOrdering} <: AbstractMixedPTLattice{O}
	δt::Float64
	Nt::Int
	δτ::Float64
	Nτ::Int
	d::Int
	ordering::O

	MixedPTLattice(δt::Real, N::Int, δτ::Real, Ni::Int, d::Int, ordering::MixedFockOrdering) = new{typeof(ordering)}(
								convert(Float64, δt), N, convert(Float64, δτ), Ni, d, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
"""
    MixedPTLattice(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, d::Int=2, ordering::MixedFockOrdering=M2M1_m1M1m2M2())

Convenience constructor for a mixed PT lattice.
"""
MixedPTLattice(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, d::Int=2, ordering::MixedFockOrdering=M2M1_m1M1m2M2()) = MixedPTLattice(
							δt, Nt, δτ, Nτ, d, ordering)
Base.similar(x::MixedPTLattice; δt::Real=x.δt, Nt::Int=x.Nt, δτ::Real=x.δτ, Nτ::Int=x.Nτ, d::Int=x.d, ordering::MixedFockOrdering=x.ordering) = MixedPTLattice(
			δt, Nt, δτ, Nτ, d, ordering)

function Base.getproperty(x::MixedPTLattice, s::Symbol)
	if s == :t
		return x.Nt * x.δt
	elseif s == :β
		return x.Nτ * x.δτ
	elseif s == :T
		return 1 / x.β
	elseif s == :ts
		return 0:x.δt:x.t
	elseif s == :τs
		return 0:x.δτ:x.β
	else
		getfield(x, s)
	end
end

Base.length(x::MixedPTLattice) = 2 * x.Nt + x.Nτ

# acending order for real branch, descending order for imag time
"""
    getindex(x::MixedPTLattice{<:M2M1_m1M1m2M2}, a::BosonicContourIndex)

Map the contour index `a` (time step `a.j`, branch `:+`, `:-` or `:τ`) to the global lattice position.
"""
function Base.getindex(x::MixedPTLattice{<:M2M1_m1M1m2M2}, a::BosonicContourIndex)
	b = branch(a)
	@boundscheck begin
		(b in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if b == :τ
			(1 <= a.j <= x.Nτ) || throw(BoundsError(1:x.Nτ, a.j))
		else
			(1 <= a.j <= x.Nt) || throw(BoundsError(1:x.Nt, a.j))
		end
	end

	if b == :+
		2*(a.j-1)+2 + x.Nτ
	elseif b == :-
		2*(a.j-1)+1 + x.Nτ
	else
		x.Nτ-a.j + 1
	end
end



# key is timestep, conj, branch, band
"""
    indexmappings(lattice::MixedPTLattice)

Return a `Dict{Tuple{Int, Symbol}, Int}` mapping `(time step, branch)` to the global index of the mixed PT lattice.
"""
function indexmappings(lattice::MixedPTLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.Nτ
		f = :τ
		r[(i, f)] = lattice[BosonicContourIndex(i, f)]
	end
	for i in 1:lattice.Nt
		for f in (:+, :-)
			r[(i, f)] = lattice[BosonicContourIndex(i, f)]
		end
	end
	return r
end

