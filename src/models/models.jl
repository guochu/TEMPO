include("def.jl")
include("unitary/unitary.jl")
include("dissipative.jl")



"""
    sysdynamics(gmps::ADT, lattice::AbstractADTLattice, model::AbstractImpurityOperator, args...; kwargs...)

Apply the system dynamics of the impurity model on a copy of the ADT (calling [`sysdynamics!`](@ref)) and return the evolved ADT.
"""
sysdynamics(gmps::ADT, lattice::AbstractADTLattice, model::AbstractImpurityOperator, args...; kwargs...) = sysdynamics!(copy(gmps), lattice, model, args...; kwargs...)
"""
    sysdynamics(lattice::AbstractADTLattice, model::AbstractImpurityOperator, args...; kwargs...)

Construct an ADT starting from the vacuum state and apply the system dynamics of the impurity model.
"""
function sysdynamics(lattice::AbstractADTLattice, model::AbstractImpurityOperator, args...; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(model))
	sysdynamics!(vacuumstate(T, lattice), lattice, model, args...; kwargs...)
end 


"""
    sysdynamics(gmps::ProcessTensor, lattice::AbstractPTLattice, model::AbstractImpurityOperator; kwargs...)

Apply the system dynamics of the impurity model on a copy of the process tensor (calling [`sysdynamics!`](@ref)) and return the evolved PT.
"""
sysdynamics(gmps::ProcessTensor, lattice::AbstractPTLattice, model::AbstractImpurityOperator; kwargs...) = sysdynamics!(copy(gmps), lattice, model; kwargs...)
"""
    sysdynamics(lattice::AbstractPTLattice, model::AbstractImpurityOperator; kwargs...)

Construct a process tensor starting from the vacuum state and apply the system dynamics of the impurity model.
"""
function sysdynamics(lattice::AbstractPTLattice, model::AbstractImpurityOperator; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(model))
	return sysdynamics!(vacuumstate(T, lattice), lattice, model; kwargs...)
end 

