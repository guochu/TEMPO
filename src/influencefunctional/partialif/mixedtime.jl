function hybriddynamics!(gmps::ADT, lattice::MixedADTLattice1Order, corr::AbstractMixedCorrelationFunction, hyb::AdditiveHyb; trunc::TruncationScheme=DefaultITruncation)
	op = hyb.op
	(lattice.d == length(op)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	ds = [lattice.d for i in 1:length(lattice)]

	for b1 in branches(lattice)
		k1 = ifelse(b1==:τ, lattice.Nτ, lattice.Nt)
		for i in 1:k1
			i′ = (b1 == :τ) ? i+1 : i
			pos1 = index(lattice, i′, branch=b1)
			pos2s = Int[]
			coefs = scalartype(lattice)[]
			for b2 in branches(lattice)
				k2 = ifelse(b2==:τ, lattice.Nτ, lattice.Nt)
				for j in 1:k2
					j′ = (b1==:τ) ? j+1 : j
					pos2 = index(lattice, j′, branch=b2)
					coef = index(corr, i, j, b1=b1, b2=b2)
					push!(pos2s, pos2)
					push!(coefs, coef)
				end
			end
			tmp = partialif_densemps(ds, pos1, pos2s, op, coefs)
			mult!(gmps, tmp, trunc=trunc)
		end
	end
	return gmps
end


