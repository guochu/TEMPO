
struct PTExpectationCache{M<:ProcessTensor, G<:Tuple, L<:AbstractPTLattice, Hl, Hr} 
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

function PTExpectationCache(lattice::AbstractPTLattice, As::Tuple, args...)
	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())

	L = num_transfers(lattice)
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

environments(lattice::ImagPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor}) = PTExpectationCache(lattice, (A, B...))
environments(lattice::MixedPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor}) = PTExpectationCache(lattice, (A, B...))
environments(lattice::RealPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor}; ρ₀::AbstractMatrix=_eye(phydim(lattice))) = PTExpectationCache(lattice, (A, B...), ρ₀)

expectationvalue(m::ContourOperator, cache::PTExpectationCache) = expectationvalue(tofockprodterm(m, cache.lattice), cache)
expectationvalue(m::AbstractFockTerm, cache::PTExpectationCache) = expectation(m, cache) / Zvalue(cache)
function expectation(m::AbstractFockTerm, cache::PTExpectationCache)
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
	# println("j=", j′, " k=", k′)
	for tj in k′:-1:j′
		right = transfer_right(right, tj, tn) 
	end	
	return contract_center(left, right)
end

pos2step(lat::ImagPTLattice, pos::Int) = imag_pos2step(pos)
pos2step(lat::RealPTLattice, pos::Int) = real_pos2step(pos)
function pos2step(lat::MixedPTLattice, pos::Int)
	if pos <= lat.Nτ
		return imag_pos2step(pos)
	else
		return real_pos2step(pos - lat.Nτ) + lat.Nτ
	end
end

imag_pos2step(pos::Int) = pos
real_pos2step(pos::Int) = div(pos+1, 2) 

