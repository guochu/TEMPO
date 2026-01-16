abstract type AbstractImpurityOperator end


struct ImpurityHamiltonian{M<:AbstractMatrix} <: AbstractImpurityOperator
	m::M
end
propagator(h::ImpurityHamiltonian, lat, b::Symbol) = _get_propagator(h.m, lat, b)
propagator(h::ImpurityHamiltonian, lat; branch::Symbol=:τ) = propagator(h, lat, branch)
phydim(h::ImpurityHamiltonian) = size(h.m, 1)
ImpurityHamiltonian(d::Int) = ImpurityHamiltonian(zeros(d, d))

TO.scalartype(::Type{ImpurityHamiltonian{M}}) where {M} = scalartype(M)



# dissipative impurity
struct ImpurityLindbladian <: AbstractImpurityOperator
	m::Array{ComplexF64, 4}
end

phydim(h::ImpurityLindbladian) = size(h.m, 1)

ImpurityLindbladian(d::Int) = ImpurityLindbladian(zeros(ComplexF64, d, d, d, d))
TO.scalartype(::Type{ImpurityLindbladian}) = ComplexF64
ImpurityLindbladian(L::LindbladOperator) = ImpurityLindbladian(L.m)
ImpurityLindbladian(H::AbstractMatrix, jumpops::Vector{<:AbstractMatrix}) = ImpurityLindbladian(lindbladoperator(H, jumpops))
