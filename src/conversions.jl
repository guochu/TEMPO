

# similarptlattice(lattice::ImagADTLattice1Order) = ImagPTLattice1Order(N=lattice.N, δτ=lattice.δτ, d=lattice.d, ordering=lattice.ordering)
# similarptlattice(lattice::RealADTLattice1Order) = RealPTLattice1Order(N=lattice.N, δt=lattice.δt, d=lattice.d, ordering=lattice.ordering)
# similarptlattice(lattice::MixedADTLattice1Order) = MixedPTLattice1Order(Nt=lattice.Nt, Nτ=lattice.Nτ, δt=lattice.δt, δτ=lattice.δτ, ordering=lattice.ordering)

# similaradtlattice(lattice::ImagPTLattice1Order) = ImagADTLattice1Order(N=lattice.N, δτ=lattice.δτ, d=lattice.d, ordering=lattice.ordering)
# similaradtlattice(lattice::RealPTLattice1Order) = RealADTLattice1Order(N=lattice.N, δt=lattice.δt, d=lattice.d, ordering=lattice.ordering)
# similaradtlattice(lattice::MixedPTLattice1Order) = MixedADTLattice1Order(Nt=lattice.Nt, Nτ=lattice.Nτ, δt=lattice.δt, δτ=lattice.δτ, ordering=lattice.ordering)


# function toadt(lattice::ImagADTLattice1Order, pt::ProcessTensor; trunc::TruncationScheme=DefaultIntegrationTruncation)
# 	lattice2 = similaradtlattice(lattice)
	
# end
