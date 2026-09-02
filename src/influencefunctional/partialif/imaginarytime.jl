"""
	hybriddynamics!(gmps, lattice, corr, hyb[, alg])

Multiply the partial influence functional of each contour lattice point in-place into `gmps` to construct the full influence functional (`ADT`). Supports `AdditiveHyb` coupling on imaginary-time (`ImagADTLattice1Order`), real-time (`RealADTLattice`) and mixed-time (`MixedADTLattice1Order`) lattices; methods for `NonAdditiveHyb`/`NonDiagonalHyb` coupling on PT lattices and for the `PartialIF`/`XTRGIF` algorithms can be found in the corresponding files.

# Arguments
- `gmps::ADT`: augmented density tensor, modified in-place.
- `lattice::AbstractADTLattice`: ADT contour lattice.
- `corr::AbstractCorrelationFunction`: bath correlation function.
- `hyb::AdditiveHyb`: additive (diagonal) system-bath coupling.
- `trunc::TruncationScheme=DefaultITruncation`: truncation scheme used after each multiplication.

# Returns
The modified `gmps`.

See also [`hybriddynamics`](@ref) and [`hybriddynamics_naive`](@ref).
"""
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