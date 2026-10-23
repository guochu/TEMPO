# TDVP-based influence functional construction (finite version).
#
# TDVPIF / TDVP2IF are influence functional construction algorithms on the same
# footing as XTRGIF and PartialIF. They view the influence functional
# as the "equilibrium state" IF = exp(H) of the influence operator H (the MPO
# form returned by `influenceoperators`), and compute it with a second-order
# TDVP imaginary-time flow
#
#     dz/dτ = H·z ,   τ : 0 → 1 ,
#
# starting from the identity influence functional z(0) = I (β = 0), so that the
# result is z(1) = e^H·z(0). Each flow step is one forward-backward TDVP sweep.
# `TDVPIF` runs the single-site sweeps (the center tensor is evolved by +δτ/2
# through Krylov exponentiation of the local effective map, factorized by QR
# (LQ), and the bond matrix by -δτ/2, with a full +δτ step at the turning site)
# on a zero-padded bond profile; `TDVP2IF` runs the two-site sweeps (the site
# pair is evolved jointly by +δτ/2 and re-split by a truncating SVD under
# `alg.trunc`, followed by the single-site backward half-step -δτ/2) on the
# unpadded guess, the bond dimension adapting during the flow.
#
# The sweeps themselves are FiniteMPSAlgorithms' TDVP (`sweep!`), driven on the
# chains' payloads. The manifold (and hence the cache) is selected by the chain
# kind:
#
# * ADT: the chain is an MPS whose physical leg is the fused (o, i) pair, so the
#   generator acts pointwise — FMA's Hadamard-product flow `dz/dτ = H ∘ z`
#   (`HadamardTDVPCache` + `HadamardTDVP` / `HadamardTDVP2`);
# * ProcessTensor: the chain is a genuine MPO, so the generator acts by left
#   multiplication — FMA's density-operator flow `dz/dτ = H·z`
#   (`TDVPCache` + `TDVP1` / `TDVP2`).
#
# `stepsize` of those algorithms is the complex time increment itself, so one
# `sweep!` applies exp(stepsize·H); the flow runs `nsteps = 1/δ` sweeps of
# `stepsize = δτ = 1/nsteps`, i.e. exp(+1·H) in total.
#
# For the single-site flow the initial identity state is zero-padded up to the
# bond dimension D and canonicalized without truncation, so that the zero-weight
# directions become orthonormal directions of the environments and the sweeps
# can populate the full bond profile min(d^j, d^{L-j}, D); truncation is applied
# only in the final canonicalization. The two-site flow needs no padding: the
# pair updates grow the bonds up to the cap of `alg.trunc` (see
# `_tdvpif_prepare!`).


# on real-time lattices `influenceoperators` returns 4 branch MPOs
# ((+,+), (+,−), (−,+), (−,−)); the total influence operator driving the
# flow is their SUM (MPS/MPO direct sum, bond dimensions add up). Indeed the
# differential IF built by `influenceoperatorstepper` is the Hadamard
# (site-wise) product e^{dt·h1}∘e^{dt·h2}∘e^{dt·h3}∘e^{dt·h4}, and in the
# element-wise algebra of ADT/PT products e^a∘e^b = e^{a+b}, so the generator
# of the full IF is h1+h2+h3+h4 (NOT the Hadamard product h1∘h2∘h3∘h4).
# NOTE: the direct-sum bond dimensions of the four branches add up, so the
# branches are summed one by one, compressing with SVD canonicalization after
# each addition to keep the bond dimension bounded. The compression uses the
# tight `Defaults.alg_trunc_sys()` (not `alg.trunc`): the compression error of H
# is exponentially amplified by the flow (IF = e^H), so a loose tolerance
# would degrade the accuracy of the influence functional.
function _tdvpif_hamiltonian(lattice, corr, hyb, alg::AbstractTDVPIFAlgorithm)
	h1, h2, h3, h4 = influenceoperators(lattice, corr, hyb, algexpan=alg.alg_expan)
	orth = Orthogonalize(SVD(), Defaults.alg_trunc_sys(); normalize=false)
	H = h1 + h2
	canonicalize!(H, alg=orth)
	H = H + h3
	canonicalize!(H, alg=orth)
	H = H + h4
	canonicalize!(H, alg=orth)
	return H
end

# ======================================================================
# flow engine (FiniteMPSAlgorithms TDVP)
# ======================================================================

# promote the flow state to the scalar type of the influence operator H when the
# latter is wider (e.g. ComplexF64 from a complex Prony exponential expansion);
# no-op otherwise. `complex` exists for both ADT and ProcessTensor.
function _tdvpif_promote_flowstate(z::Dense1DTN, H::Dense1DTN)
	T = promote_type(scalartype(z), scalartype(H))
	(T <: Complex) || return z
	(T == scalartype(z)) && return z
	return complex(z)
end

# The flow must be driven by the generator's *represented value*
# (value = scaling^L · ∏tensors), but FMA's caches contract its site tensors
# only. The Hadamard cache folds `scaling(H)^L` into its local generators, while
# the density-operator cache has no such hook (a plain `MPO` carries no scaling),
# so there the scaling is absorbed into the site tensors first.
_absorb_generator_scaling!(H::ADT) = H
_absorb_generator_scaling!(H::ProcessTensor) = _absorb_scaling!(H)

# absorb the global scaling factor of a tensor network into its site tensors
# (value = scaling^L · ∏tensors  →  value = 1^L · ∏(scaling·tensors)) and reset
# the scaling to 1, so that downstream code contracting the raw site tensors
# represents the same operator regardless of the gauge
function _absorb_scaling!(x::Dense1DTN)
	sca = scaling(x)
	(sca == 1) && return x
	for i in 1:length(x)
		x[i] = sca * x[i]
	end
	setscaling!(x, 1)
	return x
end

# the flow runs on the chains' payloads (`z` is evolved in place through the
# cache's reference to it): the ADT's fused-leg MPS takes the pointwise
# (Hadamard) generator, the ProcessTensor's MPO the left-multiplying one
_tdvpif_cache(H::ADT, z::ADT) = HadamardTDVPCache(z.parent, H.parent)
_tdvpif_cache(H::ProcessTensor, z::ProcessTensor) = TDVPCache(z.parent, MPO(H.parent))

# `stepsize` is the complex time increment applied by one `sweep!`, i.e. the
# flow integrates dz/dτ = H·z over τ : 0 → 1 with exp(stepsize·H) per sweep.
# `ishermitian=false` (Arnoldi) throughout: only the bra side of the
# environments enters conjugated, so the projected generators are not hermitian
# in general. The single-site variant (`TDVPIF`) freezes the bond profile of the
# guess; the two-site variant (`TDVP2IF`) adapts it, with `trunc` capping the
# pair re-splits.
_tdvpif_stepalg(::ADT, alg::TDVPIF, δτ::Float64) = HadamardTDVP(stepsize=δτ, verbosity=alg.verbosity)
_tdvpif_stepalg(::ProcessTensor, alg::TDVPIF, δτ::Float64) =
	TDVP1(stepsize=δτ, ishermitian=false, verbosity=alg.verbosity)
_tdvpif_stepalg(::ADT, alg::TDVP2IF, δτ::Float64) =
	HadamardTDVP2(stepsize=δτ, trunc=alg.trunc, verbosity=alg.verbosity)
_tdvpif_stepalg(::ProcessTensor, alg::TDVP2IF, δτ::Float64) =
	TDVP2(stepsize=δτ, trunc=alg.trunc, ishermitian=false, verbosity=alg.verbosity)

# preparation of the flow state z before the sweeps (after complex promotion):
# the single-site sweeps resize nothing, so z is zero-padded up to the flow bond
# dimension D — `changebond!` also right-canonicalizes, which the caches require;
# the two-site sweeps grow the bonds themselves (up to the cap of `trunc`), so no
# padding is needed — only the canonical form is imposed, since the identity
# influence functional (a plain delta tensor per site) is not isometric.
_tdvpif_prepare!(z::Dense1DTN, alg::TDVPIF) = changebond!(z, alg.trunc.D)
_tdvpif_prepare!(z::Dense1DTN, alg::TDVP2IF) =
	canonicalize!(z; alg=Orthogonalize(SVD(), NoTruncation(), normalize=false))

function _tdvpif_flow!(z::Dense1DTN, H::Dense1DTN, alg::AbstractTDVPIFAlgorithm)
	_absorb_generator_scaling!(H)
	env = _tdvpif_cache(H, z)
	nsteps = round(Int, 1 / alg.δ)
	δτ = 1 / nsteps
	stepalg = _tdvpif_stepalg(z, alg, δτ)
	for n in 1:nsteps
		sweep!(env, stepalg)
		(alg.verbosity > 1) && println("TDVPIF flow step $n/$nsteps, τ = $(n * δτ)")
	end
	(alg.verbosity > 1) && println("TDVPIF flow finished, τ = 1")
	return z
end

# ======================================================================
# ADT engine
# ======================================================================

# run the TDVP flow z(τ=1) = e^H·z(0) directly on the input state z, which may
# be the identity influence functional (β = 0) or, more generally, any MPO on
# the same lattice, e.g. the pure impurity dynamics (with Lindblad dissipation)
# obtained from `sysdynamics`; the influence operator H is thereby merged into
# the impurity dynamics in a single flow.
#
# preparation: see `_tdvpif_prepare!` — the single-site flow lifts z to the flow
# bond dimension (zero-padded, right-canonical), the two-site flow only imposes
# the canonical form (the pair updates grow the bonds up to the cap of `trunc`);
# finalization: nothing — the sweeps keep the chain canonical at the flow bond
# profile, only the Schmidt values are stale (never written by the sweeps), so
# they are marked uninitialized. The global scaling factor of z is carried
# through the flow by the `_renormalize!` bookkeeping of the payload, so the
# output value is e^H·z(0) regardless of the gauge of the input.
function _tdvpif_hybriddynamics_adt!(z::ADT, H::ADT, alg::AbstractTDVPIFAlgorithm)
	# On the imaginary axis the correlation/hybridization are real-typed, yet the
	# Prony (alg_expan) exponential expansion of a real correlation can legitimately
	# return complex exponents, which makes the influence operator H complex. The
	# flow must then run in complex arithmetic; promote the flow state accordingly.
	z = _tdvpif_promote_flowstate(z, H)
	_tdvpif_prepare!(z, alg)
	_tdvpif_flow!(z, H, alg)
	unset_svectors!(z)
	return z
end

"""
	hybriddynamics(lattice::ImagADTLattice, corr::ImagCorrelationFunction, hyb::AdditiveHyb, alg::AbstractTDVPIFAlgorithm)

Construct the influence functional on an imaginary-time ADT lattice with the `TDVPIF` (single-site) or `TDVP2IF` (two-site) algorithm: evolve the identity influence functional along the imaginary-time flow dz/dτ = H·z (H = influence operator) from τ = 0 to τ = 1 with second-order TDVP sweeps.

# Returns
The influence functional, represented as an `ADT`.

See also [`TDVPIF`](@ref), [`TDVP2IF`](@ref).
"""
function hybriddynamics(lattice::ImagADTLattice, corr::ImagCorrelationFunction, hyb::AdditiveHyb, alg::AbstractTDVPIFAlgorithm)
	T = promote_type(scalartype(corr), scalartype(hyb), Float64)
	z = vacuumstate(T, lattice)
	return hybriddynamics!(z, lattice, corr, hyb, alg)
end

"""
	hybriddynamics!(gmps::ADT, lattice::ImagADTLattice, corr::ImagCorrelationFunction, hyb::AdditiveHyb, alg::AbstractTDVPIFAlgorithm)

In-place version of the `TDVPIF`/`TDVP2IF` algorithms on ADT lattices: the influence operator H drives the TDVP flow directly on `gmps`, i.e. the flow evolves `z(τ=1) = e^H·gmps` from `z(0) = gmps` in place. This allows merging the influence functional into an arbitrary initial MPO — e.g. the pure impurity dynamics (with Lindblad dissipation) obtained from `sysdynamics` — in a single flow instead of constructing the IF separately and multiplying it in afterwards.

# Returns
The modified `gmps`.
"""
function hybriddynamics!(gmps::ADT, lattice::ImagADTLattice, corr::ImagCorrelationFunction, hyb::AdditiveHyb, alg::AbstractTDVPIFAlgorithm)
	# influence operator H in MPO (fused-leg MPS) form
	H = only(influenceoperators(lattice, corr, hyb, algexpan=alg.alg_expan))
	return _tdvpif_hybriddynamics_adt!(gmps, H, alg)
end

"""
	hybriddynamics(lattice::RealADTLattice, corr::RealCorrelationFunction, hyb::AdditiveHyb, alg::AbstractTDVPIFAlgorithm)

Construct the influence functional on a real-time ADT lattice with the `TDVPIF` (single-site) or `TDVP2IF` (two-site) algorithm: the 4 branch influence operators returned by `influenceoperators` ((+,+), (+,−), (−,+), (−,−)) are summed one by one into a single influence operator H (each partial sum compressed by SVD canonicalization with `Defaults.alg_trunc_sys()`), which then drives the same TDVP imaginary-time flow as in the imaginary-time case.

# Returns
The influence functional, represented as an `ADT`.

See also [`TDVPIF`](@ref), [`TDVP2IF`](@ref).
"""
function hybriddynamics(lattice::RealADTLattice, corr::RealCorrelationFunction, hyb::AdditiveHyb, alg::AbstractTDVPIFAlgorithm)
	T = promote_type(scalartype(corr), scalartype(hyb), Float64)
	z = vacuumstate(T, lattice)
	return hybriddynamics!(z, lattice, corr, hyb, alg)
end

"""
	hybriddynamics!(gmps::ADT, lattice::RealADTLattice, corr::RealCorrelationFunction, hyb::AdditiveHyb, alg::AbstractTDVPIFAlgorithm)

In-place version of the `TDVPIF`/`TDVP2IF` algorithms on real-time ADT lattices: the influence operator H (the sum of the 4 branch operators, see [`hybriddynamics`](@ref)) drives the TDVP flow directly on `gmps`, i.e. the flow evolves `z(τ=1) = e^H·gmps` from `z(0) = gmps` in place, merging the influence functional into an arbitrary initial MPO in a single flow.

# Returns
The modified `gmps`.
"""
function hybriddynamics!(gmps::ADT, lattice::RealADTLattice, corr::RealCorrelationFunction, hyb::AdditiveHyb, alg::AbstractTDVPIFAlgorithm)
	H = _tdvpif_hamiltonian(lattice, corr, hyb, alg)
	return _tdvpif_hybriddynamics_adt!(gmps, H, alg)
end


# ======================================================================
# PT engine
# ======================================================================

# same as `_tdvpif_hybriddynamics_adt!`, for process tensors: run the TDVP
# flow z(τ=1) = e^H·z(0) directly on the input process tensor z (the identity
# influence functional or any impurity-dynamics MPO, e.g. the output of
# `sysdynamics` with Lindblad dissipation).
function _tdvpif_hybriddynamics_pt!(z::ProcessTensor, H::ProcessTensor, alg::AbstractTDVPIFAlgorithm)
	# same complex-promotion rationale as in `_tdvpif_hybriddynamics_adt!`
	z = _tdvpif_promote_flowstate(z, H)
	_tdvpif_prepare!(z, alg)
	_tdvpif_flow!(z, H, alg)
	unset_svectors!(z)
	return z
end

"""
	hybriddynamics(lattice::ImagPTLattice, corr::ImagCorrelationFunction, hyb::GeneralHybStyle, alg::AbstractTDVPIFAlgorithm)

Construct the influence functional on an imaginary-time PT lattice with the `TDVPIF` (single-site) or `TDVP2IF` (two-site) algorithm, supporting `GeneralHybStyle` (e.g. `NonAdditiveHyb`, `NonDiagonalHyb`) coupling.

# Returns
The influence functional, represented as a `ProcessTensor`.

See also [`TDVPIF`](@ref), [`TDVP2IF`](@ref).
"""
function hybriddynamics(lattice::ImagPTLattice, corr::ImagCorrelationFunction, hyb::GeneralHybStyle, alg::AbstractTDVPIFAlgorithm)
	T = promote_type(scalartype(corr), scalartype(hyb), Float64)
	z = vacuumstate(T, lattice)
	return hybriddynamics!(z, lattice, corr, hyb, alg)
end

"""
	hybriddynamics!(gmps::ProcessTensor, lattice::ImagPTLattice, corr::ImagCorrelationFunction, hyb::GeneralHybStyle, alg::AbstractTDVPIFAlgorithm)

In-place version of the `TDVPIF`/`TDVP2IF` algorithms on PT lattices: the influence operator H drives the TDVP flow directly on `gmps`, i.e. the flow evolves `z(τ=1) = e^H·gmps` from `z(0) = gmps` in place. This allows merging the influence functional into an arbitrary initial MPO — e.g. the pure impurity dynamics (with Lindblad dissipation) obtained from `sysdynamics` — in a single flow instead of constructing the IF separately and multiplying it in afterwards.

# Returns
The modified `gmps`.
"""
function hybriddynamics!(gmps::ProcessTensor, lattice::ImagPTLattice, corr::ImagCorrelationFunction, hyb::GeneralHybStyle, alg::AbstractTDVPIFAlgorithm)
	# influence operator H in MPO form
	H = only(influenceoperators(lattice, corr, hyb, algexpan=alg.alg_expan))
	return _tdvpif_hybriddynamics_pt!(gmps, H, alg)
end

"""
	hybriddynamics(lattice::RealPTLattice, corr::RealCorrelationFunction, hyb::GeneralHybStyle, alg::AbstractTDVPIFAlgorithm)

Construct the influence functional on a real-time PT lattice with the `TDVPIF` (single-site) or `TDVP2IF` (two-site) algorithm: the 4 branch influence operators returned by `influenceoperators` ((+,+), (+,−), (−,+), (−,−)) are summed one by one into a single influence operator H (each partial sum compressed by SVD canonicalization with `Defaults.alg_trunc_sys()`), which then drives the same TDVP imaginary-time flow as in the imaginary-time case.

# Returns
The influence functional, represented as a `ProcessTensor`.

See also [`TDVPIF`](@ref), [`TDVP2IF`](@ref).
"""
function hybriddynamics(lattice::RealPTLattice, corr::RealCorrelationFunction, hyb::GeneralHybStyle, alg::AbstractTDVPIFAlgorithm)
	T = promote_type(scalartype(corr), scalartype(hyb), Float64)
	z = vacuumstate(T, lattice)
	return hybriddynamics!(z, lattice, corr, hyb, alg)
end

"""
	hybriddynamics!(gmps::ProcessTensor, lattice::RealPTLattice, corr::RealCorrelationFunction, hyb::GeneralHybStyle, alg::AbstractTDVPIFAlgorithm)

In-place version of the `TDVPIF`/`TDVP2IF` algorithms on real-time PT lattices: the influence operator H (the sum of the 4 branch operators, see [`hybriddynamics`](@ref)) drives the TDVP flow directly on `gmps`, i.e. the flow evolves `z(τ=1) = e^H·gmps` from `z(0) = gmps` in place, merging the influence functional into an arbitrary initial MPO in a single flow.

# Returns
The modified `gmps`.
"""
function hybriddynamics!(gmps::ProcessTensor, lattice::RealPTLattice, corr::RealCorrelationFunction, hyb::GeneralHybStyle, alg::AbstractTDVPIFAlgorithm)
	H = _tdvpif_hamiltonian(lattice, corr, hyb, alg)
	return _tdvpif_hybriddynamics_pt!(gmps, H, alg)
end
