abstract type AbstractImpurityOperator end


"""
    ImpurityHamiltonian{M<:AbstractMatrix} <: AbstractImpurityOperator

Hamiltonian model of the impurity system, wrapping a `d × d` matrix `m`, where `d` is the physical dimension.
"""
struct ImpurityHamiltonian{M<:AbstractMatrix} <: AbstractImpurityOperator
	m::M
end
propagator(h::ImpurityHamiltonian, lat, b::Symbol) = _get_propagator(h.m, lat, b)
propagator(h::ImpurityHamiltonian, lat; branch::Symbol=:τ) = propagator(h, lat, branch)
phydim(h::ImpurityHamiltonian) = size(h.m, 1)
"""
    ImpurityHamiltonian(d::Int)

Construct a zero-matrix Hamiltonian of physical dimension `d`.
"""
ImpurityHamiltonian(d::Int) = ImpurityHamiltonian(zeros(d, d))

TO.scalartype(::Type{ImpurityHamiltonian{M}}) where {M} = scalartype(M)



# dissipative impurity
"""
    ImpurityLindbladian <: AbstractImpurityOperator

Lindblad-type dissipative model of the impurity system, wrapping a `d × d × d × d` tensor `m`
(in superoperator form, for dissipative real-time evolution).
"""
struct ImpurityLindbladian <: AbstractImpurityOperator
	m::Array{ComplexF64, 4}
end

phydim(h::ImpurityLindbladian) = size(h.m, 1)

"""
    ImpurityLindbladian(d::Int)

Construct a zero Lindblad superoperator of physical dimension `d`.
"""
ImpurityLindbladian(d::Int) = ImpurityLindbladian(zeros(ComplexF64, d, d, d, d))
TO.scalartype(::Type{ImpurityLindbladian}) = ComplexF64
"""
    ImpurityLindbladian(L::LindbladOperator)

Construct an impurity Lindblad model from a `LindbladOperator`.
"""
ImpurityLindbladian(L::LindbladOperator) = ImpurityLindbladian(L.m)
"""
    ImpurityLindbladian(H::AbstractMatrix, jumpops::Vector{<:AbstractMatrix})

Construct a Lindblad superoperator from the Hamiltonian `H` and the jump operators `jumpops`.
"""
ImpurityLindbladian(H::AbstractMatrix, jumpops::Vector{<:AbstractMatrix}) = ImpurityLindbladian(lindbladoperator(H, jumpops))
"""
    ImpurityLindbladian(h::ImpurityHamiltonian)

Construct the corresponding Lindblad model from an [`ImpurityHamiltonian`](@ref) without jump operators (purely unitary evolution).
"""
ImpurityLindbladian(h::ImpurityHamiltonian) = ImpurityLindbladian(lindbladoperator(h.m, []))