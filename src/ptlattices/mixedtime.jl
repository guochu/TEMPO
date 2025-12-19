abstract type MixedPTLattice{O<:MixedFockOrdering} <: AbstractPTLattice{O} end
TO.scalartype(::Type{<:MixedPTLattice}) = ComplexF64
branches(::Type{<:MixedPTLattice}) = (:+, :-, :τ)

"""
	struct MixedGrassmannLattice1Order <: MixedGrassmannLattice

First order splitting of the real-time contour
"""
struct MixedPTLattice1Order{O<:MixedFockOrdering} <: MixedPTLattice{O}
	δt::Float64
	Nt::Int
	δτ::Float64
	Nτ::Int
	d::Int
	ordering::O

	MixedPTLattice1Order(δt::Real, N::Int, δτ::Real, Ni::Int, d::Int, ordering::MixedFockOrdering) = new{typeof(ordering)}(
								convert(Float64, δt), N, convert(Float64, δτ), Ni, d, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
MixedPTLattice1Order(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, d::Int=2, ordering::MixedFockOrdering=M2M1_m1M1m2M2()) = MixedPTLattice1Order(
							δt, Nt, δτ, Nτ, d, ordering)
Base.similar(x::MixedPTLattice1Order; δt::Real=x.δt, Nt::Int=x.Nt, δτ::Real=x.δτ, Nτ::Int=x.Nτ, d::Int=x.d, ordering::MixedFockOrdering=x.ordering) = MixedPTLattice1Order(
			δt, Nt, δτ, Nτ, d, ordering)


function MixedPTLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return MixedPTLattice1Order(; kwargs...)
	else
		error("Second orderr MixedGrassmannLattice not implemented")
	end
end

function Base.getproperty(x::MixedPTLattice1Order, s::Symbol)
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

Base.length(x::MixedPTLattice1Order) = 2 * x.Nt + x.Nτ

# acending order for real branch, descending order for imag time
function index(x::MixedPTLattice1Order{<:M2M1_m1M1m2M2}, i::Int; branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if branch == :τ
			(1 <= i <= x.Nτ) || throw(BoundsError(1:x.Nτ, i))
		else
			(1 <= i <= x.Nt) || throw(BoundsError(1:x.Nt, i))
		end
	end

	if branch == :+
		2*(i-1)+2 + x.Nτ
	elseif branch == :-
		2*(i-1)+1 + x.Nτ
	else
		x.Nτ-i + 1
	end
end



# key is timestep, conj, branch, band
function indexmappings(lattice::MixedPTLattice1Order)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.Nτ
		f = :τ
		r[(i, f)] = index(lattice, i, branch=f)
	end
	for i in 1:lattice.Nt
		for f in (:+, :-)
			r[(i, f)] = index(lattice, i, branch=f)
		end
	end
	return r
end

