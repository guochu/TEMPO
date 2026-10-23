"""
	hybriddynamics!(gmps, lattice, corr, hyb[, alg])

Multiply the partial influence functional of each contour lattice point in-place into `gmps` to construct the full influence functional (`ADT`). Supports `AdditiveHyb` coupling on imaginary-time (`ImagADTLattice`), real-time (`RealADTLattice`) and mixed-time (`MixedADTLattice`) lattices; methods for `NonAdditiveHyb`/`NonDiagonalHyb` coupling on PT lattices and for the `PartialIF`/`XTRGIF` algorithms can be found in the corresponding files.

# Arguments
- `gmps::ADT`: augmented density tensor, modified in-place.
- `lattice::AbstractADTLattice`: ADT contour lattice.
- `corr::AbstractCorrelationFunction`: bath correlation function.
- `hyb::AdditiveHyb`: additive (diagonal) system-bath coupling.
- `trunc::TruncationScheme=Defaults.alg_trunc_if()`: truncation scheme used after each multiplication.

# Returns
The modified `gmps`.

See also [`hybriddynamics`](@ref) and [`hybriddynamics_naive`](@ref).
"""
function hybriddynamics!(gmps::ADT, lattice::ImagADTLattice, corr::ImagCorrelationFunction, hyb::AdditiveHyb; trunc::TruncationScheme=Defaults.alg_trunc_if())
	k = lattice.N
	op = hyb.op
	(lattice.d == length(op)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	ds = [lattice.d for i in 1:length(lattice)]
	for i in 1:k
		pos1 = lattice[BosonicContourIndex(i+1, :τ)]
		pos2s = Int[]
		coefs = scalartype(lattice)[]
		for j in 1:k
			pos2 = lattice[BosonicContourIndex(j+1, :τ)]
			coef = corr[BosonicContourIndex(i, :τ), BosonicContourIndex(j, :τ)]
			push!(pos2s, pos2)
			push!(coefs, coef)
		end
		tmp = partialif_densemps(ds, pos1, pos2s, op, coefs)
		# println("bond dimension of $i-th partial IF is ", bonddims(tmp))
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end