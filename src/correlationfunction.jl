function correlationfunction(bath::AbstractBath, lattice::ImagADTLattice1Order)
    # @assert lattice.β == bath.β
    (lattice.β == bath.β) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δτ(bath, N=lattice.N, δτ=lattice.δτ)
end 
correlationfunction(bath::AbstractBath, lattice::RealADTLattice1Order) = Δt(bath, N=lattice.N, t=lattice.t) 
function correlationfunction(bath::AbstractBath, lattice::MixedADTLattice1Order)
    (lattice.β == bath.β) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δm(bath, Nτ=lattice.Nτ, t=lattice.t, Nt=lattice.Nt)
end  