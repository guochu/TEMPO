# function hybriddynamics!(gmps::ADT, lattice::RealADTLattice, corr::RealCorrelationFunction, hyb::AdditiveHyb; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
# 	for i in 1:lattice.Nt, b1 in branches(lattice), band1 in 1:lattice.bands
# 		pos1 = index(lattice, i, band=band1, branch=b1)
# 		pos2s = Int[]
# 		coefs = scalartype(lattice)[]
# 		for j in 1:lattice.Nt, b2 in branches(lattice), band2 in 1:lattice.bands
# 			pos2 = index(lattice, j, band=band2, branch=b2)
# 			coef = index(corr, i, j, b1=b1, b2=b2)
# 			push!(pos2s, pos2)
# 			push!(coefs, coef)
# 		end
# 		tmp = partialif_densemps(length(lattice), pos1, pos2s, coefs)
# 		mult!(gmps, tmp, trunc=trunc)
# 	end
# 	return gmps
# end