"""
	hybriddynamics!(gmps::ADT, lattice::MixedADTLattice1Order, corr::AbstractMixedCorrelationFunction, hyb::AdditiveHyb; trunc::TruncationScheme=DefaultITruncation)

`hybriddynamics!` method for the mixed-time contour (`MixedADTLattice1Order`, containing both imaginary-time and real-time branches).

On the mixed contour the correlation-function index `i` on every branch maps to the lattice index `i+1` (the same
`index(lattice, i+1)` convention as the imaginary-time-only lattice): the branch lattice index 1 (`τ=0`, `t=0`) is the
boundary/junction point glued by `boundarycondition!` and receives no influence-functional gates, while the lattice
indices `2..kτ` (`2..kt`) carry the correlation indices `1..Nτ` (`1..Nt`). See the `MixedADTLattice1Order` docstring
for details.
"""
function hybriddynamics!(gmps::ADT, lattice::MixedADTLattice1Order, corr::AbstractMixedCorrelationFunction, hyb::AdditiveHyb; trunc::TruncationScheme=DefaultITruncation)
	op = hyb.op
	(lattice.d == length(op)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	ds = [lattice.d for i in 1:length(lattice)]

	for b1 in branches(lattice)
		k1 = ifelse(b1==:τ, lattice.Nτ, lattice.Nt)
		for i in 1:k1
			pos1 = index(lattice, i+1, branch=b1)
			pos2s = Int[]
			coefs = scalartype(lattice)[]
			for b2 in branches(lattice)
				k2 = ifelse(b2==:τ, lattice.Nτ, lattice.Nt)
				for j in 1:k2
					pos2 = index(lattice, j+1, branch=b2)
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


