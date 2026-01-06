
function correlation(lattice::Union{ImagADTLattice1Order, MixedADTLattice1Order}, model::BosonicImpurity, op::ContourOperator, mpsI::ADT; 
					trunc::TruncationScheme=DefaultKTruncation)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice)
	Z = integrate(mpsK, mpsI)
	mpsK′ = sysdynamics(lattice, model, op, trunc=trunc)
	mpsK′ = boundarycondition!(mpsK′, lattice)
	v = integrate(mpsK′, mpsI)
	return v / Z
end
function correlation(lattice::RealADTLattice1Order, model::BosonicImpurity, op::ContourOperator, mpsI::ADT, ρ0::VecOrMat=ones(lattice.d); 
						trunc::TruncationScheme=DefaultKTruncation)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρ0)
	Z = integrate(mpsK, mpsI)
	mpsK′ = sysdynamics(lattice, model, op, trunc=trunc)
	mpsK′ = boundarycondition!(mpsK′, lattice, ρ₀=ρ0)
	v = integrate(mpsK′, mpsI)
	return v / Z	
end

function correlation(lattice::Union{ImagPTLattice1Order, MixedPTLattice1Order}, op::ContourOperator, mps::ProcessTensor)
	Z = integrate(lattice, mps)
	mps2 = apply!(op, lattice, copy(mps))
	return integrate(lattice, mps2) / Z
end

function correlation(lattice::RealPTLattice1Order, op::ContourOperator, mps::ProcessTensor, ρ0::AbstractMatrix)
	tmp = initialstate!(copy(mps), lattice, ρ0)
	Z = integrate(lattice, tmp)
	mps2 = apply!(op, lattice, copy(mps))
	mps2 = initialstate!(mps2, lattice, ρ0)
	return integrate(lattice, mps2) / Z
end
