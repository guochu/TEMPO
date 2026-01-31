

struct PTExpectationCache{M<:ProcessTensor, G<:Tuple, L<:_AllowdFinitePTLattices, Hl, Hr} 
	A::M
	Bs::G
	lattice::L
	hleft::Hl
	hright::Hr	
end

Base.length(x::PTExpectationCache) = length(x.lattice)
Zvalue(x::PTExpectationCache) = contract_center(x.hleft[1], x.hright[1])
Zvalue2(x::PTExpectationCache) = contract_center(x.hleft[end], x.hright[end])
leftenv(x::PTExpectationCache, j::Int) = x.hleft[j]
rightenv(x::PTExpectationCache, j::Int) = x.hright[j+1]

function PTExpectationCache(lattice::_AllowdFinitePTLattices, As::Tuple, args...)
	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())

	L = lattice.N
	m = PTTransferMatrix(lattice, As...)

	left = l_LL(m)
	hleft = Vector{typeof(left)}(undef, L+1)
	hleft[1] = left
	for i in 1:L
		left = transfer_left(left, i, m)
		hleft[i+1] = left
	end
	right = r_RR(m, args...)
	hright = Vector{typeof(right)}(undef, L+1)
	hright[L+1] = right
	for i in L:-1:1
		right = transfer_right(right, i, m) 
		hright[i] = right
	end
	return PTExpectationCache(first(As), Base.tail(As), lattice, hleft, hright)
end

environments(lattice::ImagPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor}) = PTExpectationCache(lattice, (A, B...))
environments(lattice::RealPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor}; ρ₀::AbstractMatrix=_eye(phydim(lattice))) = PTExpectationCache(lattice, (A, B...), ρ₀)

expectation(m::ContourOperator, cache::PTExpectationCache) = expectation(tofockprodterm(m, cache.lattice), cache)
expectation(m::AbstractFockTerm, cache::PTExpectationCache) = expectationvalue(m, cache) / Zvalue(cache)
function expectationvalue(m::AbstractFockTerm, cache::PTExpectationCache)
	# println("positions ", m.positions)
	j, k = m.positions[1], m.positions[end]
	j′ = pos2step(cache.lattice, j)
	k′ = pos2step(cache.lattice, k)
	@assert j′ <= k′
	left = leftenv(cache, j′)  
	right = rightenv(cache, k′) 
	As = (cache.A, cache.Bs...)
	A2 = apply!(m, copy(last(As)))
	As2 = Base.front(As)
	tn = PTTransferMatrix(cache.lattice, As2..., A2)
	for tj in k′:-1:j′
		right = transfer_right(right, tj, tn) 
	end	
	return contract_center(left, right)
end

pos2step(lat::ImagPTLattice, pos::Int) = pos
pos2step(lat::RealPTLattice, pos::Int) = div(pos+1, 2) 
