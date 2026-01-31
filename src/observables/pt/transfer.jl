struct PTTransferMatrix{_LAT <: AbstractPTLattice, T<:Number, N}
	lattice::_LAT
	scaling::Float64
	states::NTuple{N, Vector{Array{T, 4}}}
end

Base.length(x::PTTransferMatrix) = length(x.states[1])
TO.scalartype(::Type{PTTransferMatrix{_L, T, N}}) where {_L, T, N} = T
TEMPO.scaling(x::PTTransferMatrix) = x.scaling

function PTTransferMatrix(lattice::AbstractPTLattice, states::NTuple{N, Vector{Array{T, 4}}}; scaling::Real=1) where {T, N}
	L = length(states[1])
	for i in 2:N
		(length(states[i]) == L) || throw(DimensionMismatch("input states size mismatch"))
	end
	return PTTransferMatrix(lattice, float(scaling), states)
end 
function PTTransferMatrix(lattice::_AllowdPTLattices, states::Vararg{ProcessTensor}; scaling::Real=scaling(states...))
	(length(states[1]) == length(lattice)) || throw(DimensionMismatch("lattice size mismatch with ProcessTensor size"))
	return PTTransferMatrix(lattice, map(x->x.data, states), scaling=scaling)
end 


# the meaning for imag and real ptlattices are very different!!
# only support imag and real PT lattices currently
# support for mixed PT lattices is a bit tricky
const _AllowdFinitePTLattices{O} = Union{ImagPTLattice{O}, RealPTLattice{O}} where {O}

l_LL(f, m::PTTransferMatrix{L, T, N}) where {L<:_AllowdFinitePTLattices, T, N} = f(T, ntuple(i->space_l(m.states[i][1]), N)..., phydim(m.lattice), phydim(m.lattice))
r_RR(f, m::PTTransferMatrix{L, T, N}) where {L<:_AllowdFinitePTLattices, T, N} = f(T, ntuple(i->space_r(m.states[i][end]), N)..., phydim(m.lattice), phydim(m.lattice))


function l_LL(m::PTTransferMatrix{<:_AllowdFinitePTLattices, T, N}) where {T, N}
	d = phydim(m.lattice)
	a = _eye(T, d)
	return reshape(a, ntuple(x->1, N)..., d, d)
end 
function r_RR(m::PTTransferMatrix{<:_AllowdFinitePTLattices, T, N}) where {T, N}
	d = phydim(m.lattice)
	a = _eye(T, d)
	return reshape(a, ntuple(x->1, N)..., d, d)
end
# initial state
function r_RR(m::PTTransferMatrix{<:RealPTLattice, T, N}, ρ₀::AbstractMatrix) where {T, N}
	d = phydim(m.lattice)
	ρ = convert(Matrix{T}, ρ₀)
	return reshape(ρ, ntuple(x->1, N)..., d, d)
end

# nperiod(m::PTTransferMatrix) = div(length(m), length(m.lattice))

# random_boundaries(m::PTTransferMatrix) = random_left_boundary(m), random_right_boundary(m)
# random_left_boundary(m::PTTransferMatrix) = l_LL(randn, m)
# random_right_boundary(m::PTTransferMatrix) = r_RR(randn, m)

function transfer_left(left::AbstractArray, m::PTTransferMatrix) 
	for i in 1:m.lattice.k
		left = transfer_left(left, i, m.lattice, scaling(m), m.states...)
	end
	return left
end
function transfer_right(m::PTTransferMatrix, right::AbstractArray) 
	for i in m.lattice.k:-1:1
		right = transfer_right(right, i, m.lattice, scaling(m), m.states...)
	end
	return right
end

# imaginary time 
function transfer_left(left::DenseMPSTensor, i::Int, lattice::ImagPTLattice, sca, x::Vector{<:DenseMPOTensor})
	@tensor tmp[3,4,5] := left[1,2,5] * x[i][1,2,3,4]
	return lmul!(sca, tmp)
end
function transfer_left(left::DenseMPOTensor, i::Int, lattice::ImagPTLattice, sca, x::Vector{<:DenseMPOTensor}, y::Vector{<:DenseMPOTensor})
	@tensor tmp[4,6,7,8] := left[1,2,3,8] * x[i][1,3,4,5] * y[i][2,5,6,7]
	return lmul!(sca, tmp)
end

function transfer_right(right::DenseMPOTensor, i::Int, lattice::ImagPTLattice, sca, x::Vector{<:DenseMPOTensor})
	@tensor tmp[1,2,5,6] := x[i][1,2,3,4] * right[3,4,5,6]
	return lmul!(sca, tmp)
end
function transfer_right(right::DenseMPOTensor, i::Int, lattice::ImagPTLattice, sca, x::Vector{<:DenseMPOTensor}, y::Vector{<:DenseMPOTensor})
	@tensor tmp[1,7,2] := x[i][1,2,3,4] * right[3,5,6] * y[i][7,4,5,6]
	return lmul!(sca, tmp)
end

# real time
function transfer_left(left::DenseMPSTensor, i::Int, lattice::RealPTLattice, sca, x::Vector{<:DenseMPOTensor})
	@tensor tmp[6,5,7] := left[1,2, 3] * x[2i-1][1,2,4,5] * x[2i][4,3,6,7]
	return lmul!(sca^2, tmp)
end
function transfer_left(left::DenseMPOTensor, i::Int, lattice::RealPTLattice, sca, x::Vector{<:DenseMPOTensor}, y::Vector{<:DenseMPOTensor})
	@tensor tmp[5,7,4,8] := left[1,2, 3,4] * x[2i-1][1,3,5,6] * y[2i-1][2,6,7,8]
	@tensor left[5,7,4,8] := tmp[1,2, 3,4] * x[2i][1,3,5,6] * y[2i][2,6,7,8]
	return lmul!(sca^2, left)
end
function transfer_right(right::DenseMPSTensor, i::Int, lattice::RealPTLattice, sca, x::Vector{<:DenseMPOTensor})
	@tensor tmp[6,7,5] := right[1,2,3] * x[2i][4,5,1,3] * x[2i-1][6,7,4,2]
	return lmul!(sca^2, tmp)
end
function transfer_right(right::DenseMPOTensor, i::Int, lattice::RealPTLattice, sca, x::Vector{<:DenseMPOTensor}, y::Vector{<:DenseMPOTensor})
	@tensor tmp[7,5,8,3] := right[1,2,3,4] * y[2i][5,6,2,4] * x[2i][7,8,1,6]
	@tensor right[7,5,8,3] := tmp[1,2,3,4] * y[2i-1][5,6,2,4] * x[2i-1][7,8,1,6]
	return lmul!(sca^2, right)
end

scaling(x::ProcessTensor, y::ProcessTensor, zs::ProcessTensor...) = scaling(x) * scaling(y) * prod(map(scaling, zs))
