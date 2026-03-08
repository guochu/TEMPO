
# struct PTMixedExpectationCache{M<:ProcessTensor, G<:Tuple, L<:MixedPTLattice, Hl, Hr} 
# 	A::M
# 	Bs::G
# 	lattice::L
# 	hleft::Hl
# 	hright::Hr	
# end

# Base.length(x::PTMixedExpectationCache) = length(x.lattice)
# Zvalue(x::PTMixedExpectationCache) = contract_center(x.hleft[1], x.hright[1])
# Zvalue2(x::PTMixedExpectationCache) = contract_center(x.hleft[end], x.hright[end])
# leftenv(x::PTMixedExpectationCache, j::Int) = x.hleft[j]
# rightenv(x::PTMixedExpectationCache, j::Int) = x.hright[j+1]

# function PTMixedExpectationCache(lattice::MixedPTLattice, As::Tuple)
# 	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())

# 	L = lattice.Nτ + lattice.Nt
# 	m = PTTransferMatrix(lattice, As...)

# 	left = l_LL(m)
# 	hleft = Vector{typeof(left)}(undef, L+1)
# 	hleft[1] = left
# 	for i in 1:lattice.Nτ
# 		left = transfer_imag_left(left, i, m)
# 		hleft[i+1] = left
# 	end
# 	for i in 1:lattice.Nt
# 		left = transfer_real_left(left, i, m)
# 		hleft[lattice.Nτ+1+i] = left
# 	end
# 	right = r_RR(m)
# 	hright = Vector{typeof(right)}(undef, L+1)
# 	hright[L+1] = right
# 	for i in lattice.Nt:-1:1
# 		right = transfer_real_right(right, i, m) 
# 		hright[lattice.Nτ + i] = right
# 	end
# 	for i in lattice.Nτ:-1:1
# 		right = transfer_imag_right(right, i, m)
# 		hright[i] = right
# 	end
# 	return PTMixedExpectationCache(first(As), Base.tail(As), lattice, hleft, hright)
# end

# environments(lattice::MixedPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor}) = PTMixedExpectationCache(lattice, (A, B...))

# expectationvalue(m::ContourOperator, cache::PTMixedExpectationCache) = expectationvalue(tofockprodterm(m, cache.lattice), cache)
# expectationvalue(m::AbstractFockTerm, cache::PTMixedExpectationCache) = expectation(m, cache) / Zvalue(cache)
# function expectation(m::AbstractFockTerm, cache::PTMixedExpectationCache)
# 	# println("positions ", m.positions)
# 	j, k = m.positions[1], m.positions[end]
# 	j′ = pos2step(cache.lattice, j)
# 	k′ = pos2step(cache.lattice, k)
# 	@assert j′ <= k′
# 	left = leftenv(cache, j′)  
# 	right = rightenv(cache, k′) 
# 	As = (cache.A, cache.Bs...)
# 	A2 = apply!(m, copy(last(As)))
# 	As2 = Base.front(As)
# 	tn = PTTransferMatrix(cache.lattice, As2..., A2)
# 	for tj in k′:-1:j′
# 		right = transfer_right(right, tj, tn) 
# 	end	
# 	return contract_center(left, right)
# end

# function pos2step(lat::MixedPTLattice, pos::Int)
# 	if pos <= lattice.Nτ
# 		return imag_pos2step(pos)
# 	else
# 		return real_pos2step(pos - lattice.Nτ) + lattice.Nτ
# 	end
# end




# function transfer_left(left::AbstractArray, m::PTTransferMatrix{<:MixedPTLattice}) 
# 	left = transfer_imag_left(left, m)
# 	return transfer_real_left(left, m)
# end
# function transfer_right(m::PTTransferMatrix{<:MixedPTLattice}, right::AbstractArray) 
# 	right = transfer_real_right(m, right)
# 	return transfer_imag_right(m, right)
# end

# function transfer_real_right(m::PTTransferMatrix{<:MixedPTLattice}, right::AbstractArray) 
# 	for i in m.lattice.Nτ:-1:1
# 		right = transfer_real_right(right, i, m)
# 	end
# 	return right
# end


# function transfer_real_left(left::AbstractArray, m::PTTransferMatrix{<:MixedPTLattice}) 
# 	for i in 1:m.lattice.Nτ
# 		left = transfer_real_left(left, i, m)
# 	end
# 	return left
# end
# function transfer_real_right(m::PTTransferMatrix{<:MixedPTLattice}, right::AbstractArray) 
# 	for i in m.lattice.Nτ:-1:1
# 		right = transfer_real_right(right, i, m)
# 	end
# 	return right
# end

# function transfer_imag_left(left::AbstractArray, m::PTTransferMatrix{<:MixedPTLattice}) 
# 	for i in 1:m.lattice.Nτ
# 		left = transfer_imag_left(left, i, m)
# 	end
# 	return left
# end
# function transfer_imag_right(m::PTTransferMatrix{<:MixedPTLattice}, right::AbstractArray) 
# 	for i in m.lattice.Nτ:-1:1
# 		right = transfer_imag_right(right, i, m)
# 	end
# 	return right
# end

function transfer_left(left::AbstractArray, i::Int, m::PTTransferMatrix{<:MixedPTLattice}) 
	if i <= m.lattice.Nτ
		return transfer_imag_left(left, i, m)
	else
		return transfer_real_left(left, i-m.lattice.Nτ, m)
	end
end
function transfer_right(right::AbstractArray, i::Int, m::PTTransferMatrix{<:MixedPTLattice})
	if i <= m.lattice.Nτ
		return transfer_imag_right(right, i, m)
	else
		return transfer_real_right(right, i-m.lattice.Nτ, m)
	end
end

transfer_imag_left(left::AbstractArray, i::Int, m::PTTransferMatrix{<:MixedPTLattice}) = transfer_imag_left(left, i, m.lattice, scaling(m), m.states...)
transfer_imag_right(right::AbstractArray, i::Int, m::PTTransferMatrix{<:MixedPTLattice}) = transfer_imag_right(right, i, m.lattice, scaling(m), m.states...)

transfer_real_left(left::AbstractArray, i::Int, m::PTTransferMatrix{<:MixedPTLattice}) = transfer_real_left(left, i, m.lattice, scaling(m), m.states...)
transfer_real_right(right::AbstractArray, i::Int, m::PTTransferMatrix{<:MixedPTLattice}) = transfer_real_right(right, i, m.lattice, scaling(m), m.states...)


# mixed time
function transfer_imag_left(left::DenseMPSTensor, i::Int, lattice::MixedPTLattice, sca, x::Vector{<:DenseMPOTensor})
	@assert i <= lattice.Nτ
	@tensor tmp[4,5,3] := left[1,2,3] * x[i][1,2,4,5]
	return lmul!(sca, tmp)
end
function transfer_imag_left(left::DenseMPOTensor, i::Int, lattice::MixedPTLattice, sca, x::Vector{<:DenseMPOTensor}, y::Vector{<:DenseMPOTensor})
	@assert i <= lattice.Nτ
	@tensor tmp[5,7,8,4] := left[1,2,3,4] * x[i][1,3,5,6] * y[i][2,6,7,8]
	return lmul!(sca, tmp)
end
function transfer_real_left(left::DenseMPSTensor, i::Int, lattice::MixedPTLattice, sca, x::Vector{<:DenseMPOTensor})
	@assert i <= lattice.Nt
	@tensor tmp[7,4,6] := left[1,2, 3] * x[2i-1+lattice.Nτ][1,4,5,2] * x[2i+lattice.Nτ][5,6,7,3]
	return lmul!(sca^2, tmp)
end
function transfer_real_left(left::DenseMPOTensor, i::Int, lattice::MixedPTLattice, sca, x::Vector{<:DenseMPOTensor}, y::Vector{<:DenseMPOTensor})
	@assert i <= lattice.Nt
	@tensor tmp[6,8,4,5] := left[1,2,3,4] * x[2i-1+lattice.Nτ][1,5,6,7] * y[2i-1+lattice.Nτ][2,7,8,3]
	@tensor left[6,8,4,5] := tmp[1,2,3,4] * x[2i+lattice.Nτ][1,5,6,7] * y[2i+lattice.Nτ][2,7,8,3]
	return lmul!(sca^2, left)
end


function transfer_imag_right(right::DenseMPSTensor, i::Int, lattice::MixedPTLattice, sca, x::Vector{<:DenseMPOTensor})
	@assert i <= lattice.Nτ
	@tensor tmp[1,2,5] := x[i][1,2,3,4] * right[3,4,5]
	return lmul!(sca, tmp)
end
function transfer_imag_right(right::DenseMPOTensor, i::Int, lattice::MixedPTLattice, sca, x::Vector{<:DenseMPOTensor}, y::Vector{<:DenseMPOTensor})
	@assert i <= lattice.Nτ
	@tensor tmp[1,8,2,7] := x[i][1,2,3,4] * right[3,5,6,7] * y[i][8,4,5,6]
	return lmul!(sca, tmp)
end
function transfer_real_right(right::DenseMPSTensor, i::Int, lattice::MixedPTLattice, sca, x::Vector{<:DenseMPOTensor})
	@assert i <= lattice.Nt
	@tensor tmp[6,7,5] := right[1,2,3] * x[2i+lattice.Nτ][4,3,1,5] * x[2i-1+lattice.Nτ][6,2,4,7]
	return lmul!(sca^2, tmp)
end
function transfer_real_right(right::DenseMPOTensor, i::Int, lattice::MixedPTLattice, sca, x::Vector{<:DenseMPOTensor}, y::Vector{<:DenseMPOTensor})
	@assert i <= lattice.Nt
	@tensor tmp[8,5,7,3] := right[1,2,3,4] * y[2i+lattice.Nτ][5,6,2,7] * x[2i+lattice.Nτ][8,4,1,6] 
	@tensor right[8,5,7,3] := tmp[1,2,3,4] * y[2i-1+lattice.Nτ][5,6,2,7] * x[2i-1+lattice.Nτ][8,4,1,6] 
	return lmul!(sca^2, right)
end

