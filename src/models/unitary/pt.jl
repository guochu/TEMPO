# PT

# process tensor (nonadditive baths)
"""
    sysdynamics!(gmps::ProcessTensor, lattice::ImagPTLattice, model::ImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation)

In-place version of [`sysdynamics`](@ref): apply the bare impurity Hamiltonian `model` to the given process tensor, evolving along the imaginary-time branch.
"""
function sysdynamics!(gmps::ProcessTensor, lattice::ImagPTLattice, model::ImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation)
	return sysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
end 


"""
    sysdynamics!(gmps::ProcessTensor, lattice::RealPTLattice, model::ImpurityHamiltonian;
                 branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)

Apply system dynamics in place to a process tensor on the real-time contour: with `branch=nothing`, acts successively on the `:+` and `:-` branches;
when `branch` is specified, acts only on the corresponding branch.
"""
function sysdynamics!(gmps::ProcessTensor, lattice::RealPTLattice, model::ImpurityHamiltonian; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		return sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? sysdynamics_forward!(gmps, lattice, model; trunc=trunc) : sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	end
end 

"""
    sysdynamics!(gmps::ProcessTensor, lattice::MixedPTLattice, model::ImpurityHamiltonian;
                 branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)

Apply system dynamics in place to a process tensor on the mixed contour: with `branch=nothing`, acts successively on the forward, backward and imaginary-time branches;
when `branch` (`:+`, `:-` or `:τ`) is specified, acts only on the corresponding branch.
"""
function sysdynamics!(gmps::ProcessTensor, lattice::MixedPTLattice, model::ImpurityHamiltonian; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
		return sysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
	else
		if branch == :+
			return sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		elseif branch == :-
			return sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
		else
			(branch == :τ) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
			return sysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
		end
	end
end 


sysdynamics_forward!(mps::ProcessTensor, lattice::AbstractPTLattice, model::ImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :+, lattice.Nt; trunc=trunc)
sysdynamics_backward!(mps::ProcessTensor, lattice::AbstractPTLattice, model::ImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :-, lattice.Nt; trunc=trunc)
sysdynamics_imaginary!(mps::ProcessTensor, lattice::AbstractPTLattice, model::ImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :τ, lattice.Nτ; trunc=trunc)


function _sysdynamics_util!(gmps::ProcessTensor, lattice::AbstractPTLattice, model::ImpurityHamiltonian, branch::Symbol, N::Int; trunc::TruncationScheme=DefaultKTruncation)
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