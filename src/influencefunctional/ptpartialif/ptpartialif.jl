include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")

hybriddynamics(gmps::ProcessTensor, lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, bs::NonAdditiveHyb; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr, bs; kwargs...)
function hybriddynamics(lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, bs::NonAdditiveHyb; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(bs), scalartype(corr))
	return hybriddynamics!(vacuumstate(T, lattice), lattice, corr, bs; kwargs...)
end 


hybriddynamics_naive(gmps::ProcessTensor, lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, bs::NonAdditiveHyb; kwargs...) = hybriddynamics_naive!(copy(gmps), lattice, corr, bs; kwargs...)
function hybriddynamics_naive(lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, bs::NonAdditiveHyb; kwargs...) 
	T = promote_type(scalartype(lattice), scalartype(bs), scalartype(corr))
	hybriddynamics_naive!(vacuumstate(T, lattice), lattice, corr, bs; kwargs...)
end


# naive implementation with N^2 gate operations
# function hybriddynamics_naive!(gmps::ProcessTensor, lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, hyb::NonAdditiveHyb; trunc::TruncationScheme=DefaultITruncation)
# 	z = hyb.op
# 	(lattice.d == size(z, 1) == size(z, 2)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
# 	d = lattice.d
# 	z2 = z * z
# 	zz = kron(z, z)
# 	orth = Orthogonalize(SVD(), trunc)

# 	for b1 in branches(lattice)
# 		k1 = (b1 == :τ) ? lattice.Nτ : lattice.Nt
# 		for i in 1:k1
# 			tmp = vacuumstate(lattice)
# 			ind1 = ContourIndex(i, b1) 
# 			for b2 in branches(lattice)
# 				k2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
# 				for j in 1:k2
# 					coef = index(corr, i, j, b1=b1, b2=b2)
# 					ind2 = ContourIndex(j, b2) 
# 					if ind1 == ind2
# 						m = exp(coef .* z2)
# 						# t = ContourOperator(ind1, m)
# 						t = FockTermS(lattice[ind1], m)
# 					else
# 						m = exp(coef .* zz)
# 						# t = ContourOperator([ind1, ind2], [z, z])
# 						t = FockTermS((lattice[ind1], lattice[ind2]), reshape(m, (d,d,d,d)))
# 					end
# 					apply!(t, tmp)
# 					canonicalize!(tmp, alg=orth)
# 				end
# 			end
# 			gmps = mult!(gmps, tmp, trunc=trunc)			
# 		end
# 	end
# 	return gmps
# end

function hybriddynamics_naive!(gmps::ProcessTensor, lattice::AbstractPTLattice, corr::AbstractCorrelationFunction, hyb::NonAdditiveHyb; 
								trunc::TruncationScheme=DefaultITruncation)
	for b1 in branches(lattice)
		k1 = (b1 == :τ) ? lattice.Nτ : lattice.Nt
		for i in 1:k1
			ind1 = ContourIndex(i, branch=b1)
			tmp = partialif_naive(lattice, ind1, corr, hyb, trunc=trunc)
			gmps = mult!(gmps, tmp, trunc=trunc)			
		end
	end
	return gmps
end