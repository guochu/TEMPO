function hybriddynamics!(gmps::ProcessTensor, lattice::ImagPTLattice1Order, corr::ImagCorrelationFunction, hyb::NonAdditiveHyb; trunc::TruncationScheme=DefaultITruncation)
	k = lattice.N
	op = hyb.op
	(lattice.d == size(op, 1)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	orth = Orthogonalize(SVD(), trunc)
	for i in 1:k
		pos1 = index(lattice, i)
		pos2s = Int[]
		coefs = scalartype(lattice)[]
		for j in 1:k
			pos2 = index(lattice, j)
			coef = index(corr, i, j)
			push!(pos2s, pos2)
			push!(coefs, coef)
		end
		tmp = partialif_densempo(pos1, pos2s, op, coefs)
		# println("bond dimension of $i-th partial IF is ", bond_dimensions(tmp))
		apply!(tmp, gmps)
		canonicalize!(gmps, alg=orth)
	end
	return gmps
end