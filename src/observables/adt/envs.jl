
struct ADTExpectationCache{L<:AbstractADTLattice, M<:AbstractImpurityOperator, _I<:ADT, _Hl, _Hr} 
	lattice::L
	model::M
	mpsI::_I
	hleft::_Hl
	hright::_Hr
end

Base.length(x::ADTExpectationCache) = length(x.lattice)
TO.scalartype(::Type{ADTExpectationCache{L, M, _I, _Hl, _Hr}}) where {L, M, _I, _Hl, _Hr} = promote_type(scalartype(L), scalartype(M), scalartype(_I))
Zvalue(x::ADTExpectationCache) = TO.scalar(x.hleft[end])
leftenv(x::ADTExpectationCache, j::Int) = x.hleft[j]
rightenv(x::ADTExpectationCache, j::Int) = x.hright[j+1]

function environments(lattice::AbstractADTLattice, model::AbstractImpurityOperator, mpsI::ADT; trunc::TruncationScheme=DefaultKTruncation)
	(length(lattice) == lattice(mpsI)) || throw(DimensionMismatch("lattice size and IF size mismatch"))
	(phydim(lattice) == phydim(model)) || throw(DimensionMismatch("lattice phydim mismatch with model phydim"))
	L = length(lattice)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	left = l_LL(mpsK, mpsI)
	hleft = Vector{typeof(left)}(undef, L+1)
	hleft[1] = left
	for i in 1:L
		left = left * TransferMatrix(i, mpsK, mpsI)
		hleft[i+1] = left
	end
	right = r_RR(mpsK, mpsI)
	hright = Vector{typeof(right)}(undef, L+1)
	hright[L+1] = right
	for i in L:-1:1
		right = TransferMatrix(i, mpsK, mpsI) * right
		hright[i] = right
	end
	return ADTExpectationCache(lattice, model, mpsI, hleft, hright)
end