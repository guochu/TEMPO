"""
	hybriddynamics!(gmps::ADT, lattice::RealADTLattice, corr::RealCorrelationFunction, hyb::AdditiveHyb; trunc::TruncationScheme=Defaults.alg_trunc_if())

`hybriddynamics!` method for the real-time contour (`RealADTLattice`).
"""
function hybriddynamics!(gmps::ADT, lattice::RealADTLattice, corr::RealCorrelationFunction, hyb::AdditiveHyb; trunc::TruncationScheme=Defaults.alg_trunc_if())
	op = hyb.op
	(lattice.d == length(op)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	ds = [lattice.d for i in 1:length(lattice)]

	for i in 1:lattice.Nt, b1 in branches(lattice)
		pos1 = lattice[BosonicContourIndex(i, b1)]
		pos2s = Int[]
		coefs = scalartype(lattice)[]
		for j in 1:lattice.Nt, b2 in branches(lattice)
			pos2 = lattice[BosonicContourIndex(j, b2)]
			coef = corr[BosonicContourIndex(i, b1), BosonicContourIndex(j, b2)]
			push!(pos2s, pos2)
			push!(coefs, coef)
		end
		tmp = partialif_densemps(ds, pos1, pos2s, op, coefs)
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end