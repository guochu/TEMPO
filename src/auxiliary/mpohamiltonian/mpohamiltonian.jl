"""
	struct MPOHamiltonian{M <: AbstractSparseMPOTensor}

A generic MPO which stores a chain of AbstractSparseMPOTensor (Matrix of MPOTensors)

For finite system, the first site tensor is understood as the first row of the 
first AbstractSparseMPOTensor, and the last site tensor is understood as the last 
column of the last AbstractSparseMPOTensor
"""
struct MPOHamiltonian{M <: AbstractSparseMPOTensor}
	data::Vector{M}

function MPOHamiltonian{M}(data::AbstractVector) where {M <: AbstractSparseMPOTensor}
	@assert !isempty(data) 
	(size(data[1], 1) == size(data[end], 2) ) || throw(DimensionMismatch())
	for i in 1:length(data)-1
		(size(data[i], 2) == size(data[i+1], 1)) || throw(DimensionMismatch())

	end
	new{M}(convert(Vector{M}, data))
end

end

Base.length(x::MPOHamiltonian) = length(x.data)
Base.getindex(x::MPOHamiltonian, i::Int) = getindex(x.data, i)
Base.firstindex(x::MPOHamiltonian) = firstindex(x.data)
Base.lastindex(x::MPOHamiltonian) = lastindex(x.data)
Base.copy(x::MPOHamiltonian) = MPOHamiltonian(copy(x.data))

TO.scalartype(::Type{MPOHamiltonian{M}}) where {M} = scalartype(M)

MPOHamiltonian(data::AbstractVector{M}) where {M <: AbstractSparseMPOTensor} = MPOHamiltonian{M}(data)
MPOHamiltonian(data::Vector{Matrix{Any}}) = MPOHamiltonian([SparseMPOTensor(item) for item in data])


function Base.getindex(x::MPOHamiltonian, i::Int, j::Int, k::Int)
	x[i][j, k]
end 

Base.getindex(x::MPOHamiltonian, i::Colon, j::Int, k::Int) = [getindex(x, i,j,k) for i in 1:length(x)]


# """
# 	MPO(h::MPOHamiltonian, L::Int) 
	
# Conversion of an MPOHamiltonian into a finite dense MPO
# """
# MPO(h::MPOHamiltonian{<:SchurMPOTensor}) = _tompo(h, 1, size(h[end], 2))
# MPO(h::MPOHamiltonian{<:SparseMPOTensor}; rowl::Int=1, colr::Int=1) = _tompo(h, rowl, colr)

function _tompo(h::MPOHamiltonian, leftrow::Int, rightcol::Int) 
	L = length(h)
	(L >= 2) || throw(ArgumentError("size of MPO must at least be 2"))
	# isstrict(h) || throw(ArgumentError("only strict MPOHamiltonian is allowed"))
	T = scalartype(h)

	mpotensors = Vector{Array(T, 4)}(undef, L)
	dj = phydim(h[1])

	tmp = zeros(T, 1, dj, size(h[1], 2), dj)
	for i in 1:size(h1[1], 2)
		tmp[1, :, i, :] = h[1, leftrow, i]
	end
	mpotensors[1] = tmp
	for n in 2:L-1
		dj = phydim(h[n])
		sl, sr = size(h[n])
		tmp = zeros(T, sl, dj, sr, dj)
		for i in 1:sl, j in 1:sr
			tmp[i, :, j, :] = h[n, i, j]
		end
		mpotensors[n] = tmp
	end
	dj = phydim(h[L])
	sl = size(h[L], 1)
	tmp = zeros(T, sl, dj, 1, dj)
	# _a = size(h[L], 2)
	for i in 1:sl
		tmp[i, :, 1, :] = h[L, i, rightcol]
	end
	mpotensors[L] = tmp
	return mpotensors
end
