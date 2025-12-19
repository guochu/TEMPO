function hybriddynamics!(gmps::ADT, lattice::ImagADTLattice1Order, corr::ImagCorrelationFunction, hyb::AdditiveHyb; trunc::TruncationScheme=DefaultITruncation)
	k = lattice.N
	op = hyb.op
	(lattice.d == length(op)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	ds = [lattice.d for i in 1:length(lattice)]
	for i in 1:k
		pos1 = index(lattice, i+1)
		pos2s = Int[]
		coefs = scalartype(lattice)[]
		for j in 1:k
			pos2 = index(lattice, j+1)
			coef = index(corr, i, j)
			push!(pos2s, pos2)
			push!(coefs, coef)
		end
		tmp = partialif_densemps(ds, pos1, pos2s, op, coefs)
		# println("bond dimension of $i-th partial IF is ", bond_dimensions(tmp))
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end