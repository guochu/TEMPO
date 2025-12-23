function correlationfunction(bath::AbstractBath, lattice::Union{ImagADTLattice1Order, ImagPTLattice1Order})
    # @assert lattice.β == bath.β
    (lattice.β == bath.β) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δτ(bath, N=lattice.N, δτ=lattice.δτ)
end 
correlationfunction(bath::AbstractBath, lattice::Union{RealADTLattice1Order, RealPTLattice1Order}) = Δt(bath, N=lattice.N, t=lattice.t) 
function correlationfunction(bath::AbstractBath, lattice::Union{MixedADTLattice1Order, MixedPTLattice1Order})
    (lattice.β == bath.β) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δm(bath, Nτ=lattice.Nτ, t=lattice.t, Nt=lattice.Nt)
end  

TO.scalartype(::Type{<:ImagCorrelationFunction{<:AbstractMatrix{T}}}) where {T} = T
TO.scalartype(::Type{<:RealCorrelationFunction}) = ComplexF64
TO.scalartype(::Type{<:MixedCorrelationFunction}) = ComplexF64