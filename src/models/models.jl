abstract type AbstractBosonicImpurityHamiltonian end

# interfaces
sysdynamics(gmps::ADT, lattice::AbstractADTLattice, model::AbstractBosonicImpurityHamiltonian, args...; kwargs...) = sysdynamics!(copy(gmps), lattice, model, args...; kwargs...)
function sysdynamics(lattice::AbstractADTLattice, model::AbstractBosonicImpurityHamiltonian, args...; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(model))
	sysdynamics!(vacuumstate(T, lattice), lattice, model, args...; kwargs...)
end 


sysdynamics(gmps::ProcessTensor, lattice::AbstractPTLattice, model::AbstractBosonicImpurityHamiltonian; kwargs...) = sysdynamics!(copy(gmps), lattice, model; kwargs...)
function sysdynamics(lattice::AbstractPTLattice, model::AbstractBosonicImpurityHamiltonian; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(model))
	return sysdynamics!(vacuumstate(T, lattice), lattice, model; kwargs...)
end 


"""
	sysdynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractBosonicImpurityHamiltonian; kwargs...)

The inplace version of the function sysdynamics, which applys the bare impurity onto the given GMPS
"""
function sysdynamics!(gmps::ADT, lattice::ImagADTLattice, model::AbstractBosonicImpurityHamiltonian, args...; trunc::TruncationScheme=DefaultKTruncation)
	return sysdynamics_imaginary!(gmps, lattice, model, args...; trunc=trunc)
end 


function sysdynamics!(gmps::ADT, lattice::RealADTLattice, model::AbstractBosonicImpurityHamiltonian, args...; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward!(gmps, lattice, model, args...; trunc=trunc)
		return sysdynamics_backward!(gmps, lattice, model, args...; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? sysdynamics_forward!(gmps, lattice, model, args...; trunc=trunc) : sysdynamics_backward!(gmps, lattice, model, args...; trunc=trunc)
	end
end 

function sysdynamics!(gmps::ADT, lattice::MixedADTLattice, model::AbstractBosonicImpurityHamiltonian, args...; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward!(gmps, lattice, model, args...; trunc=trunc)
		sysdynamics_backward!(gmps, lattice, model, args...; trunc=trunc)
		return sysdynamics_imaginary!(gmps, lattice, model, args...; trunc=trunc)
	else
		if branch == :+
			return sysdynamics_forward!(gmps, lattice, model, args...; trunc=trunc)
		elseif branch == :-
			return sysdynamics_backward!(gmps, lattice, model, args...; trunc=trunc)
		else
			(branch == :τ) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
			return sysdynamics_imaginary!(gmps, lattice, model, args...; trunc=trunc)
		end
	end
end 



sysdynamics_forward!(gmps::ADT, lattice::AbstractADTLattice, model::AbstractBosonicImpurityHamiltonian, args...; kwargs...) = error("sysdynamics_forward! not implemented for model $(typeof(model))")
sysdynamics_backward!(gmps::ADT, lattice::AbstractADTLattice, model::AbstractBosonicImpurityHamiltonian, args...; kwargs...) = error("sysdynamics_backward! not implemented for model $(typeof(model))")
sysdynamics_imaginary!(gmps::ADT, lattice::AbstractADTLattice, model::AbstractBosonicImpurityHamiltonian, args...; kwargs...) = error("sysdynamics_imaginary! not implemented for model $(typeof(model))")



# process tensor (nonadditive baths)
function sysdynamics!(gmps::ProcessTensor, lattice::ImagPTLattice, model::AbstractBosonicImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation)
	return sysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
end 


function sysdynamics!(gmps::ProcessTensor, lattice::RealPTLattice, model::AbstractBosonicImpurityHamiltonian; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		return sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? sysdynamics_forward!(gmps, lattice, model; trunc=trunc) : sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	end
end 

function sysdynamics!(gmps::ProcessTensor, lattice::MixedPTLattice, model::AbstractBosonicImpurityHamiltonian; 
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



sysdynamics_forward!(gmps::ProcessTensor, lattice::AbstractPTLattice, model::AbstractBosonicImpurityHamiltonian; kwargs...) = error("sysdynamics_forward! not implemented for model $(typeof(model))")
sysdynamics_backward!(gmps::ProcessTensor, lattice::AbstractPTLattice, model::AbstractBosonicImpurityHamiltonian; kwargs...) = error("sysdynamics_backward! not implemented for model $(typeof(model))")
sysdynamics_imaginary!(gmps::ProcessTensor, lattice::AbstractPTLattice, model::AbstractBosonicImpurityHamiltonian; kwargs...) = error("sysdynamics_imaginary! not implemented for model $(typeof(model))")


include("sysdynamics.jl")
include("offdiagobs.jl")