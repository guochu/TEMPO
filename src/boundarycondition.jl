boundarycondition(x::FockMPS, lattice::AbstractFockLattice; kwargs...) = boundarycondition!(copy(x), lattice; kwargs...)


function boundarycondition!(x::FockMPS, lattice::ImagFockLattice1Order; trunc::TruncationScheme=DefaultIntegrationTruncation)
	pos1, pos2 = index(lattice, 1), index(lattice, lattice.k)
	d = lattice.d
	apply!(FockTerm((pos1, pos2), _eye(d)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end


# function boundarycondition!(x::FockMPS, lattice::RealFockLattice1Order; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
# 	d = lattice.phydims[band]
# 	pos1, pos2 = index(lattice, 1, band=band, branch=:+), index(lattice, lattice.N, band=band, branch=:+)
# 	apply!(FockTerm((pos1, pos2), _eye(d)), x)
# 	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
# 	pos1, pos2 = index(lattice, lattice.N, band=band, branch=:-), index(lattice, 1, band=band, branch=:-)
# 	apply!(FockTerm((pos1, pos2), _eye(d)), x)
# 	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
# 	return x
# end


# function boundarycondition!(x::FockMPS, lattice::MixedFockLattice1Order; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
# 	pos1, pos2 = index(lattice, lattice.Nt, band=band, branch=:-), index(lattice, lattice.Nt, band=band, branch=:+)
# 	apply!(FockTerm((pos1, pos2), _eye(d)), x)

# 	pos1, pos2 = index(lattice, 1, band=band, branch=:τ), index(lattice, 1, band=band, branch=:-)
# 	apply!(FockTerm((pos1, pos2), _eye(d)), x)
	
# 	pos1, pos2 = index(lattice, 1, band=band), index(lattice, lattice.Nτ, band=band, branch=:τ)
# 	apply!(FockTerm((pos1, pos2), _eye(d)), x)
# 	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
# 	return x
# end