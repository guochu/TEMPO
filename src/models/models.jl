include("def.jl")
include("unitary/unitary.jl")
include("dissipative.jl")



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

