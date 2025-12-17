abstract type ImagFockLattice{O<:ImagFockOrdering} <: AbstractFockLattice{O} end
TO.scalartype(::Type{<:ImagFockLattice}) = Float64
branches(::Type{<:ImagFockLattice}) = (:τ,)
TimeOrderingStyle(x::ImagFockLattice) = ImaginaryTimeOrderingStyle(x)

struct ImagFockLattice1Order{O<:ImagFockOrdering} <: ImagFockLattice{O}
	δτ::Float64
	d::Int
	N::Int
	ordering::O

	ImagFockLattice1Order(δτ::Real, d::Int, N::Int, ordering::ImagFockOrdering) = new{typeof(ordering)}(float(δτ), d, N, ordering)
end

ImagFockLattice1Order(; δτ::Real, N::Int, d::Int=2, ordering::ImagFockOrdering=M1N1()) = ImagFockLattice1Order(δτ, d, N, ordering)
Base.similar(x::ImagFockLattice1Order; δτ::Real=x.δτ, d::Int=x.d, N::Int=x.N, ordering::ImagFockOrdering=x.ordering) = ImagFockLattice1Order(δτ, d, N, ordering)
# similargrassmannlattice(x::ImagFockLattice1Order, δτ::Real=x.δτ, bands::Int=x.bands, N::Int=x.N, 
# 						ordering::ImagGrassmannOrdering=similargrassmannordering(x.ordering)) = GrassmannLattice(contour=:imag, δτ=δτ, N=N, bands=bands, ordering=ordering)

function ImagFockLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return ImagFockLattice1Order(; kwargs...)
	else
		error("Second orderr ImagGrassmannLattice not implemented")
	end
end

Base.length(x::ImagFockLattice) = x.k
function Base.getproperty(x::ImagFockLattice1Order, s::Symbol)
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


function index(x::ImagFockLattice{<:M2M1}, i::Int; branch::Symbol=:τ)
	@boundscheck begin
		(1 <= i <= x.k) || throw(BoundsError(1:x.k, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	return x.k-i+1
end



function indexmappings(lattice::ImagFockLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.k
		r[(i, :τ)] = index(lattice, i)
	end
	return r
end
