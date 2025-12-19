abstract type ImagPTLattice{O<:ImagFockOrdering} <: AbstractPTLattice{O} end
TO.scalartype(::Type{<:ImagPTLattice}) = Float64
branches(::Type{<:ImagPTLattice}) = (:τ,)
TimeOrderingStyle(x::ImagPTLattice) = ImaginaryTimeOrderingStyle(x)

struct ImagPTLattice1Order{O<:ImagFockOrdering} <: ImagPTLattice{O}
	δτ::Float64
	d::Int
	N::Int
	ordering::O

	ImagPTLattice1Order(δτ::Real, d::Int, N::Int, ordering::ImagFockOrdering) = new{typeof(ordering)}(float(δτ), d, N, ordering)
end

ImagPTLattice1Order(; δτ::Real, N::Int, d::Int=2, ordering::ImagFockOrdering=M2M1()) = ImagPTLattice1Order(δτ, d, N, ordering)
Base.similar(x::ImagPTLattice1Order; δτ::Real=x.δτ, d::Int=x.d, N::Int=x.N, ordering::ImagFockOrdering=x.ordering) = ImagPTLattice1Order(δτ, d, N, ordering)

function ImagPTLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return ImagPTLattice1Order(; kwargs...)
	else
		error("Second orderr ImagGrassmannLattice not implemented")
	end
end

Base.length(x::ImagPTLattice) = x.N
function Base.getproperty(x::ImagPTLattice1Order, s::Symbol)
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


function index(x::ImagPTLattice{<:M2M1}, i::Int; branch::Symbol=:τ)
	@boundscheck begin
		(1 <= i <= x.N) || throw(BoundsError(1:x.N, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	return x.N-i+1
end



function indexmappings(lattice::ImagPTLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.N
		r[(i, :τ)] = index(lattice, i)
	end
	return r
end
