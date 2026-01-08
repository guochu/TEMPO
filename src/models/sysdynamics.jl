abstract type AbstractImpurityOperator end

# interfaces
sysdynamics(gmps::ADT, lattice::AbstractADTLattice, model::AbstractImpurityOperator, args...; kwargs...) = sysdynamics!(copy(gmps), lattice, model, args...; kwargs...)
function sysdynamics(lattice::AbstractADTLattice, model::AbstractImpurityOperator, args...; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(model))
	sysdynamics!(vacuumstate(T, lattice), lattice, model, args...; kwargs...)
end 


sysdynamics(gmps::ProcessTensor, lattice::AbstractPTLattice, model::AbstractImpurityOperator; kwargs...) = sysdynamics!(copy(gmps), lattice, model; kwargs...)
function sysdynamics(lattice::AbstractPTLattice, model::AbstractImpurityOperator; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(model))
	return sysdynamics!(vacuumstate(T, lattice), lattice, model; kwargs...)
end 





# sysdynamics_forward!(gmps::ADT, lattice::AbstractADTLattice, model::AbstractImpurityOperator, args...; kwargs...) = error("sysdynamics_forward! not implemented for model $(typeof(model))")
# sysdynamics_backward!(gmps::ADT, lattice::AbstractADTLattice, model::AbstractImpurityOperator, args...; kwargs...) = error("sysdynamics_backward! not implemented for model $(typeof(model))")
# sysdynamics_imaginary!(gmps::ADT, lattice::AbstractADTLattice, model::AbstractImpurityOperator, args...; kwargs...) = error("sysdynamics_imaginary! not implemented for model $(typeof(model))")






# sysdynamics_forward!(gmps::ProcessTensor, lattice::AbstractPTLattice, model::AbstractImpurityOperator; kwargs...) = error("sysdynamics_forward! not implemented for model $(typeof(model))")
# sysdynamics_backward!(gmps::ProcessTensor, lattice::AbstractPTLattice, model::AbstractImpurityOperator; kwargs...) = error("sysdynamics_backward! not implemented for model $(typeof(model))")
# sysdynamics_imaginary!(gmps::ProcessTensor, lattice::AbstractPTLattice, model::AbstractImpurityOperator; kwargs...) = error("sysdynamics_imaginary! not implemented for model $(typeof(model))")


include("unitary.jl")
include("dissipative.jl")