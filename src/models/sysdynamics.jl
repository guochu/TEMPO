struct BosonicImpurity{M<:AbstractMatrix} <: AbstractBosonicImpurityHamiltonian
	m::M
end
propagator(h::BosonicImpurity, lat, b::Symbol) = _get_propagator(h.m, lat, b)
propagator(h::BosonicImpurity, lat; branch::Symbol=:τ) = propagator(h, lat, branch)
phydim(h::BosonicImpurity) = size(h.m, 1)

TO.scalartype(::Type{BosonicImpurity{M}}) where {M} = scalartype(M)

# # Ĥ = Ωσ̂ₓ
# spinboson(;Ω::Real=0) = BosonicImpurity(Ω .* pauli_x())


sysdynamics_forward!(mps::ADT, lattice::AbstractADTLattice, model::BosonicImpurity, args...; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :+, lattice.Nt, args...; trunc=trunc)
sysdynamics_backward!(mps::ADT, lattice::AbstractADTLattice, model::BosonicImpurity, args...; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :-, lattice.Nt, args...; trunc=trunc)
sysdynamics_imaginary!(mps::ADT, lattice::AbstractADTLattice, model::BosonicImpurity, args...; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :τ, lattice.Nτ, args...; trunc=trunc)

function _sysdynamics_util!(gmps::ADT, lattice::AbstractADTLattice, model::BosonicImpurity, branch::Symbol, N::Int; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	U = propagator(model, lattice, branch)
	# data = decompose_to_mps(U)
	alg = Orthogonalize(SVD(), trunc)
	for j in 1:N
		a, b = (branch == :-) ? (j, j+1) : (j+1, j)
        pos1, pos2 = index(lattice, a, branch=branch), index(lattice, b, branch=branch)
        t = ADTTerm((pos1, pos2), U)
        apply!(t, gmps)
        canonicalize!(gmps, alg=alg)			
	end
	return gmps
end


function _get_propagator(h, lattice, b::Symbol)
	if b == :τ
		return exp(-lattice.δτ .* h)
	elseif b == :+
		return exp(-im*lattice.δt .* h)
	else
		return exp(im*lattice.δt .* h)
	end
end


sysdynamics_forward!(mps::ProcessTensor, lattice::AbstractPTLattice, model::BosonicImpurity, args...; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :+, lattice.Nt, args...; trunc=trunc)
sysdynamics_backward!(mps::ProcessTensor, lattice::AbstractPTLattice, model::BosonicImpurity, args...; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :-, lattice.Nt, args...; trunc=trunc)
sysdynamics_imaginary!(mps::ProcessTensor, lattice::AbstractPTLattice, model::BosonicImpurity, args...; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :τ, lattice.Nτ, args...; trunc=trunc)


function _sysdynamics_util!(gmps::ProcessTensor, lattice::AbstractPTLattice, model::BosonicImpurity, branch::Symbol, N::Int; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	U = propagator(model, lattice, branch)
	# data = decompose_to_mps(U)
	alg = Orthogonalize(SVD(), trunc)
	for j in 1:N
        t = ContourOperator(ContourIndex(j, branch), U)
        apply!(t, lattice, gmps)			
	end
	canonicalize!(gmps, alg=alg)
	return gmps
end


# function pauli_matrices()
# 	s_SP = Array{Float64, 2}([0 0; 1 0])
# 	s_SM = Array{Float64, 2}([0 1; 0 0])
# 	s_Z = Array{Float64, 2}([-1 0; 0 1])
# 	s_x = s_SP+s_SM
# 	s_y = -im*(s_SP-s_SM)
# 	return s_x, s_y, s_Z
# end

pauli_x() = Matrix{Float64}([0 1; 1 0])
pauli_y() = Matrix{ComplexF64}([0 im; -im 0])
pauli_z() = Matrix{Float64}([-1 0; 0 1])
spin_up() = Matrix{Float64}([0 0; 0 1])
spin_down() = Matrix{Float64}([1 0; 0 0])