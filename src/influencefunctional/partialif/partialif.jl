# abstract type PartialIFStyle end

# struct TriangularPartialIF <: PartialIFStyle end
# struct RectangularPartialIF <: PartialIFStyle end


include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


hybriddynamics(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr, bs; kwargs...)
hybriddynamics(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...) = hybriddynamics!(vacuumstate(lattice), lattice, corr, bs; kwargs...)


hybriddynamics_naive(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...) = hybriddynamics_naive!(copy(gmps), lattice, corr, bs; kwargs...)
hybriddynamics_naive(lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, bs::AdditiveHyb; kwargs...) = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr, bs; kwargs...)



# naive implementation with N^2 gate operations
function hybriddynamics_naive!(gmps::ADT, lattice::AbstractADTLattice, corr::AbstractCorrelationFunction, hyb::AdditiveHyb; trunc::TruncationScheme=DefaultITruncation)
	z = hyb.op
	(lattice.d == length(z)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	d = lattice.d
	z2 = z .* z
	zz = reshape(kron(z, z), d, d)
	orth = Orthogonalize(SVD(), trunc)

	for b1 in branches(lattice)
		k1 = (b1 == :τ) ? lattice.Nτ : lattice.Nt
		for i in 1:k1
			tmp = vacuumstate(lattice)
			i′ = (b1 == :τ) ? i+1 : i
			pos1 = index(lattice, i′, branch=b1)
			for b2 in branches(lattice)
				k2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
				for j in 1:k2
					coef = index(corr, i, j, b1=b1, b2=b2)
					j′ = (b1==:τ) ? j+1 : j
					pos2 = index(lattice, j′, branch=b2)
					if pos1 == pos2
						m = exp.(coef .* z2)
						t = ADTTerm((pos1, ), (m, ))
					else
						m = exp.(coef .* zz)
						t = ADTTerm((pos1, pos2), m)
					end
					apply!(t, tmp)
					canonicalize!(tmp, alg=orth)
				end
			end
			gmps = mult!(gmps, tmp, trunc=trunc)			
		end
	end
	return gmps
end