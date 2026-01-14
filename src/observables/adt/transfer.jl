
struct ADTTransferMatrix{T<:Number, N}
	scaling::Float64
	states::NTuple{N, Vector{Array{T, 3}}}
end

function ADTTransferMatrix(states::NTuple{N, Vector{Array{T, 3}}}, scaling::Real) where {N, T}
	(N > 0) || throw(ArgumentError("no element"))
	L = length(states[1])
	all(x->length(x)==L, states) || throw(ArgumentError("elements have different lengths"))
	return ADTTransferMatrix(convert(Float64, scaling), states)
end



Base.length(x::ADTTransferMatrix) = length(x.states[1])
TO.scalartype(::Type{ADTTransferMatrix{T, N}}) where {T, N} = T
scaling(x::ADTTransferMatrix) = x.scaling

function transfer_left end
function transfer_right end


function transfer_left(left::Vector{<:Number}, j::Int, x::Vector{<:DenseMPSTensor})
	@tensor tmp[2,3] := left[1] * x[j][1,2,3]
	return dropdims(sum(tmp, dims=1), dims=1)
end

function transfer_right(right::Vector{<:Number}, j::Int, x::Vector{<:DenseMPSTensor})
	@tensor tmp[1,2] := x[j][1,2,3] * right[3]
	return dropdims(sum(tmp, dims=2), dims=2)
end

function transfer_left(left::Matrix{<:Number}, j::Int, x::Vector{<:DenseMPSTensor}, y::Vector{<:DenseMPSTensor})
	@tensor tmp[5,4] := left[1,2] * y[j][2,3,4] * x[j][1,3,5]
	return tmp
end

function transfer_right(right::Matrix{<:Number}, j::Int, x::Vector{<:DenseMPSTensor}, y::Vector{<:DenseMPSTensor})
	@tensor tmp[1,5] := x[j][1,2,3] * right[3,4] * y[j][5,2,4]
	return tmp
end

function Base.:*(left::Array{<:Number, N}, m::ADTTransferMatrix{<:Number, N}) where {N}
	for i in 1:length(m)
		left = lmul!(scaling(m), transfer_left(left, i, m.states...)) 
	end
	return left
end
function Base.:*(m::ADTTransferMatrix{<:Number, N}, right::Array{<:Number, N}) where {N}
	for i in length(m):-1:1
		right = lmul!(scaling(m), transfer_right(right, i, m.states...)) 
	end
	return right
end

l_LL(f, m::ADTTransferMatrix{T, N}) where {T, N} = f(T, ntuple(i->space_l(m.states[i][1]), N))
r_RR(f, m::ADTTransferMatrix{T, N}) where {T, N} = f(T, ntuple(i->space_r(m.states[i][end]), N))
l_LL(m::ADTTransferMatrix{T, N}) where {T, N} = l_LL(ones, m)
r_RR(m::ADTTransferMatrix{T, N}) where {T, N} = r_RR(ones, m)

TransferMatrix(states::Vararg{M, N}) where {M <: ADT, N} = ADTTransferMatrix(map(x->x.data, states), scaling(states...))
TransferMatrix(j::Int, states::Vararg{M, N}) where {M <: ADT, N} = ADTTransferMatrix(map(x->[x[j]], states), scaling(states...))


scaling(x::ADT, y::ADT, zs::ADT...) = scaling(x) * scaling(y) * prod(map(scaling, zs))