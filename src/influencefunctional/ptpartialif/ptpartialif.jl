include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")



# naive implementation with N^2 gate operations
function hybriddynamics_naive!(gmps::ProcessTensor, lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, hyb::NonAdditiveHyb; trunc::TruncationScheme=DefaultITruncation)
	z = hyb.op
	(lattice.d == size(z, 1) == size(z, 2)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	d = lattice.d
	z2 = z * z
	orth = Orthogonalize(SVD(), trunc)

	for b1 in branches(lattice)
		k1 = (b1 == :τ) ? lattice.Nτ : lattice.Nt
		for i in 1:k1
			tmp = vacuumstate(lattice)
			ind1 = ContourIndex(i, b1) 
			for b2 in branches(lattice)
				k2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
				for j in 1:k2
					coef = index(corr, i, j, b1=b1, b2=b2)
					ind2 = ContourIndex(j, b2) 
					if ind1 == ind2
						m = exp.(coef .* z2)
						t = ContourOperator(ind1, m)
					else
						m = exp.(coef .* zz)
						t = ContourOperator([ind1, ind2], [z, z])
					end
					apply!(t, lattice, tmp)
					canonicalize!(tmp, alg=orth)
				end
			end
			gmps = mult!(gmps, tmp, trunc=trunc)			
		end
	end
	return gmps
end