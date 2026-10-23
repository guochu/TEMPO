# Cross-check of the dense influence functionals built by QuAPI.Dense against
# TEMPO's own (untruncated) ADT / PT influence functionals, on the three
# contours. The dense tensors must agree element by element.

using QuAPI.Dense
using QuAPI: BosonicContourIndex
using FiniteMPSAlgorithms: NoTruncation

# dense coefficient tensor of an ADT (MPS): axes = lattice sites 1..L
function _tempo_adt_array(x::TEMPO.ADT)
    L = length(x)
    d = size(x[1], 2)
    return reshape(TEMPO.todense(x.parent; order=:lsb), ntuple(_ -> d, L)...)
end

# dense coefficient tensor of a ProcessTensor (MPO):
# axes = (ket_1, bra_1, ket_2, bra_2, ...), lattice site order
function _tempo_pt_array(x::TEMPO.ProcessTensor)
    L = length(x)
    d = size(x[1], 2)
    M = TEMPO.todense(x.parent; order=:lsb)
    T = reshape(M, ntuple(_ -> d, L)..., ntuple(_ -> d, L)...)
    return permutedims(T, vcat([[p, L + p] for p in 1:L]...))
end

# QuAPI.Dense basis corresponding to the PT site order
function _pt_basis(lat)
    L = length(lat)
    if lat isa ImagPTLattice
        return [BosonicContourIndex(lat.N - p + 1, :τ) for p in 1:L]
    elseif lat isa RealPTLattice
        return [BosonicContourIndex(lat.N - div(p - 1, 2), isodd(p) ? :+ : :-)
                for p in 1:L]
    else
        b = [BosonicContourIndex(lat.Nτ - p + 1, :τ) for p in 1:lat.Nτ]
        for r in 1:(L - lat.Nτ)
            push!(b, BosonicContourIndex(div(r + 1, 2), isodd(r) ? :- : :+))
        end
        return b
    end
end

@testset "QuAPI dense IF cross-check" begin
    β = 1.0
    spec = DiscreteSpectrum([1.0, 2.0], [0.25, 0.09])
    bath = bosonicbath(spec, β=β)
    λ = [1.0, -1.0]
    op = [1.0 0.3; 0.3 -1.0]
    rtol = 1.0e-11

    # ---------------- imaginary time ----------------
    N, δτ = 3, β / 3
    lat = ADTLattice(N=N, δτ=δτ, contour=:imag, d=2)
    corr = correlationfunction(bath, lat)
    IF = hybriddynamics(lat, corr, AdditiveHyb(λ); trunc=NoTruncation())
    corrq = bosonic_Δτ(spec; β=β, N=N)
    @test Array(selectdim(_tempo_adt_array(IF), N + 1, 1)) ≈
          reshape(BosonicIF(BosonicHybridizationMatrix(corrq), λ).data, ntuple(_ -> 2, N)...) rtol = rtol

    latP = PTLattice(N=N, δτ=δτ, contour=:imag, d=2)
    IFP = hybriddynamics(latP, correlationfunction(bath, latP), NonAdditiveHyb(op);
                         trunc=NoTruncation())
    AP = _tempo_pt_array(IFP)
    @test AP ≈ reshape(BosonicIF(BosonicHybridizationMatrix(corrq; split=true), op).data, ntuple(_ -> 2, 2N)...) rtol = rtol
    QP = BosonicPT(_pt_basis(latP), 2, vec(AP))
    @test QuAPI.Dense.Zvalue(QP) ≈ integrate(latP, IFP) rtol = 1.0e-10

    # ---------------- real time ----------------
    Nr, δt = 3, 0.1
    corrq = bosonic_Δt(spec; β=β, N=Nr - 1, t=(Nr - 1) * δt)   # N gated points
    latr = ADTLattice(N=Nr, δt=δt, contour=:real, d=2)
    IFr = hybriddynamics(latr, correlationfunction(bath, latr), AdditiveHyb(λ);
                         trunc=NoTruncation())
    Ar = _tempo_adt_array(IFr)
    @test Array(Ar[1, 1, :, :, :, :, :, :]) ≈
          reshape(BosonicIF(BosonicHybridizationMatrix(corrq), λ).data, ntuple(_ -> 2, 2Nr)...) rtol = rtol

    latrP = PTLattice(N=Nr, δt=δt, contour=:real, d=2)
    IFrP = hybriddynamics(latrP, correlationfunction(bath, latrP), NonAdditiveHyb(op);
                          trunc=NoTruncation())
    ArP = _tempo_pt_array(IFrP)
    @test ArP ≈ reshape(BosonicIF(BosonicHybridizationMatrix(corrq; split=true), op).data, ntuple(_ -> 2, 4Nr)...) rtol = rtol
    QPr = BosonicPT(_pt_basis(latrP), 2, vec(ArP))
    @test QuAPI.Dense.Zvalue(QPr) ≈ integrate(latrP, IFrP) rtol = 1.0e-10

    # ---------------- mixed time ----------------
    Nt, Nτ, δtm = 2, 3, 0.1
    δτm = β / Nτ
    corrq = bosonic_Δm(spec; β=β, Nτ=Nτ, t=(Nt - 1) * δtm, Nt=Nt - 1)
    latm = ADTLattice(Nt=Nt, δt=δtm, Nτ=Nτ, δτ=δτm, contour=:mixed, d=2)
    IFm = hybriddynamics(latm, correlationfunction(bath, latm), AdditiveHyb(λ);
                         trunc=NoTruncation())
    Am = _tempo_adt_array(IFm)
    Amg = Am[:, :, :, 1, 1, 1, :, :, :, :]
    # TEMPO's real-branch order is (i,-), (i,+); QuAPI's is (i,+), (i,-)
    peri = (1, 2, 3, 5, 4, 7, 6)
    @test Array(Amg) ≈
          permutedims(reshape(BosonicIF(BosonicHybridizationMatrix(corrq), λ).data,
                              ntuple(_ -> 2, Nτ + 2Nt)...), peri) rtol = rtol

    latmP = PTLattice(Nt=Nt, δt=δtm, Nτ=Nτ, δτ=δτm, contour=:mixed, d=2)
    IFmP = hybriddynamics(latmP, correlationfunction(bath, latmP), NonAdditiveHyb(op);
                          trunc=NoTruncation())
    AmP = _tempo_pt_array(IFmP)
    ax = vcat([[2p - 1, 2p] for p in peri]...)
    @test AmP ≈
          permutedims(reshape(BosonicIF(BosonicHybridizationMatrix(corrq; split=true), op).data,
                              ntuple(_ -> 2, 2 * (Nτ + 2Nt))...), ax) rtol = rtol
    QPm = BosonicPT(_pt_basis(latmP), 2, vec(AmP))
    @test QuAPI.Dense.Zvalue(QPm) ≈ integrate(latmP, IFmP) rtol = 1.0e-10
end
