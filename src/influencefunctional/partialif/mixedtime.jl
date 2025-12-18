# function hybriddynamics!(gmps::ADT, lattice::MixedADTLattice1Order, corr::AbstractMixedCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
# 	for b1 in branches(lattice)
# 		k1 = ifelse(b1==:τ, lattice.Nτ, lattice.Nt)
# 		for i in 1:k1, band1 in 1:lattice.bands
# 			pos1 = index(lattice, i, band=band1, branch=b1)
# 			pos2s = Int[]
# 			coefs = scalartype(lattice)[]
# 			for b2 in branches(lattice)
# 				k2 = ifelse(b2==:τ, lattice.Nτ, lattice.Nt)
# 				for j in 1:k2, band2 in 1:lattice.bands
# 					pos2 = index(lattice, j, band=band2, branch=b2)
# 					coef = index(corr, i, j, b1=b1, b2=b2)
# 					push!(pos2s, pos2)
# 					push!(coefs, coef)
# 				end
# 			end
# 			tmp = partialif_densemps(length(lattice), pos1, pos2s, coefs)
# 			mult!(gmps, tmp, trunc=trunc)
# 		end
# 	end
# 	return gmps
# end


