

struct ADTExpectationCache{M<:ADT, G<:Tuple, L<:AbstractADTLattice, Hl, Hr} 
	A::M
	Bs::G
	lattice::L
	hleft::Hl
	hright::Hr	
end

Base.length(x::ADTExpectationCache) = length(x.lattice)
Zvalue(x::ADTExpectationCache) = only(x.hleft[end])
leftenv(x::ADTExpectationCache, j::Int) = x.hleft[j]
rightenv(x::ADTExpectationCache, j::Int) = x.hright[j+1]

function ADTExpectationCache(lattice::AbstractADTLattice, As::Tuple)
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
	return ADTExpectationCache(first(As), Base.tail(As), lattice, hleft, hright)
end

environments(lattice::AbstractADTLattice, A::ADT, B::Vararg{ADT}) = ADTExpectationCache(lattice, (A, B...))

expectation(m::ADTTerm, cache::ADTExpectationCache) = expectationvalue(m, cache) / Zvalue(cache)
function expectationvalue(m::ADTTerm, cache::ADTExpectationCache)
	j, k = m.positions[1], m.positions[end]
	left = leftenv(cache, j)  
	right = rightenv(cache, k) 
	A2 = apply!(m, copy(cache.A))
	for tj in k:-1:j
		right = TransferMatrix(tj, A2, cache.Bs...) * right 
	end	
	# @tensor r = left[1] * right[1]
	# return r
	return contract_center(left, right)
end

function contract_center(left::Vector, right::Vector)
	@tensor r = left[1] * right[1]
	return r
end

function contract_center(left::Matrix, right::Matrix)
	@tensor r = left[1, 2] * right[1, 2]
	return r
end