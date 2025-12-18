abstract type MixedADTLattice{O<:MixedFockOrdering} <: AbstractADTLattice{O} end
TO.scalartype(::Type{<:MixedADTLattice}) = ComplexF64
branches(::Type{<:MixedADTLattice}) = (:+, :-, :τ)

"""
	struct MixedGrassmannLattice1Order <: MixedGrassmannLattice

First order splitting of the real-time contour
"""
struct MixedADTLattice1Order{O<:MixedFockOrdering} <: MixedADTLattice{O}
	δt::Float64
	Nt::Int
	δτ::Float64
	Nτ::Int
	d::Int
	ordering::O

	MixedADTLattice1Order(δt::Real, N::Int, δτ::Real, Ni::Int, d::Int, ordering::MixedFockOrdering) = new{typeof(ordering)}(
								convert(Float64, δt), N, convert(Float64, δτ), Ni, d, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
MixedADTLattice1Order(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, d::Int=2, ordering::MixedFockOrdering=M1N1_m1M1n1N1m2M2n2N2()) = MixedADTLattice1Order(
							δt, Nt, δτ, Nτ, d, ordering)
Base.similar(x::MixedADTLattice1Order; δt::Real=x.δt, Nt::Int=x.Nt, δτ::Real=x.δτ, Nτ::Int=x.Nτ, d::Int=x.d, ordering::MixedFockOrdering=x.ordering) = MixedADTLattice1Order(
			δt, Nt, δτ, Nτ, d, ordering)


function MixedADTLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return MixedADTLattice1Order(; kwargs...)
	else
		error("Second orderr MixedGrassmannLattice not implemented")
	end
end

function Base.getproperty(x::MixedADTLattice1Order, s::Symbol)
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
	elseif s == :kt
		return x.Nt + 1
	elseif s == :kτ
		return x.Nτ + 1
	else
		getfield(x, s)
	end
end

Base.length(x::MixedADTLattice1Order) = 2 * x.kt + x.kτ

# acending order for real branch, descending order for imag time
function index(x::MixedADTLattice1Order{<:M2M1_m1M1m2M2}, i::Int; branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if branch == :τ
			(1 <= i <= x.kτ) || throw(BoundsError(1:x.kτ, i))
		else
			(1 <= i <= x.kt) || throw(BoundsError(1:x.kt, i))
		end
	end

	if branch == :+
		2*(i-1)+2 + x.kτ
	elseif branch == :-
		2*(i-1)+1 + x.kτ
	else
		x.kτ-i + 1
	end
end



# key is timestep, conj, branch, band
function indexmappings(lattice::MixedADTLattice1Order)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.Nτ
		f = :τ
		r[(i, f)] = index(lattice, i, branch=f)
	end
	for i in 1:lattice.kt
		for f in (:+, :-)
			r[(i, f)] = index(lattice, i, branch=f)
		end
	end
	return r
end

