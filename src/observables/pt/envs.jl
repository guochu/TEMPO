

struct PTExpectationCache{M<:ProcessTensor, G<:Tuple, L<:_AllowdPTLattices, Hl, Hr} 
	A::M
	Bs::G
	lattice::L
	hleft::Hl
	hright::Hr	
end

Base.length(x::PTExpectationCache) = length(x.lattice)
Zvalue(x::PTExpectationCache) = only(x.hleft[end])
leftenv(x::PTExpectationCache, j::Int) = x.hleft[j]
rightenv(x::PTExpectationCache, j::Int) = x.hright[j+1]

function PTExpectationCache(lattice::_AllowdPTLattices, As::Tuple)
	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())

	L = length(lattice)

	left = l_LL(As...)
	hleft = Vector{typeof(left)}(undef, L+1)
	hleft[1] = left
	for i in 1:L
		left = left * TransferMatrix(i, As...)
		hleft[i+1] = left
	end
	right = r_RR(As...)
	hright = Vector{typeof(right)}(undef, L+1)
	hright[L+1] = right
	for i in L:-1:1
		right = TransferMatrix(i, As...) * right
		hright[i] = right
	end
	return PTExpectationCache(first(As), Base.tail(As), lattice, hleft, hright)
end

environments(lattice::ImagPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor}) = PTExpectationCache(lattice, (A, B...))
environments(lattice::RealPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor}; ρ₀::AbstractMatrix) = PTExpectationCache(lattice, (A, B...))

expectation(m::AbstractFockTerm, cache::PTExpectationCache) = expectationvalue(m, cache) / Zvalue(cache)
function expectationvalue(m::AbstractFockTerm, cache::PTExpectationCache)
	j, k = m.positions[1], m.positions[end]
	left = leftenv(cache, j)  
	right = rightenv(cache, k) 
	A2 = apply!(m, copy(cache.A))
	for tj in k:-1:j
		right = transfer_right(tj, A2, cache.Bs...) * right 
	end	
	# @tensor r = left[1] * right[1]
	# return r
	return contract_center(left, right)
end
