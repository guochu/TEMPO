function hybriddynamics!(gmps::ADT, lattice::RealADTLattice, corr::RealCorrelationFunction, hyb::AdditiveHyb; trunc::TruncationScheme=DefaultITruncation)
	op = hyb.op
	(lattice.d == length(op)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	ds = [lattice.d for i in 1:length(lattice)]

	for i in 1:lattice.Nt, b1 in branches(lattice)
		pos1 = index(lattice, i, branch=b1)
		pos2s = Int[]
		coefs = scalartype(lattice)[]
		for j in 1:lattice.Nt, b2 in branches(lattice)
			pos2 = index(lattice, j, branch=b2)
			coef = index(corr, i, j, b1=b1, b2=b2)
			push!(pos2s, pos2)
			push!(coefs, coef)
		end
		tmp = partialif_densemps(ds, pos1, pos2s, op, coefs)
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end