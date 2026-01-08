
struct ADTExpectationCache{L<:AbstractADTLattice, M<:AbstractImpurityOperator, _I<:ADT}
	lattice::L
	model::M
	mpsI::_I
end