# TDVP-based influence functional construction (finite version).
#
# TDVPIF is an influence functional construction algorithm on the same footing
# as XTRGIF and PartialIF. It views the influence functional
# as the "equilibrium state" IF = exp(H) of the influence operator H (the MPO
# form returned by `influenceoperator`), and computes it with a second-order
# single-site TDVP imaginary-time flow
#
#     dz/dτ = H·z ,   τ : 0 → 1 ,
#
# starting from the identity influence functional z(0) = I (β = 0). Each flow
# step is one forward-backward TDVP sweep: in a left-to-right (right-to-left)
# sweep the center tensor AC is evolved by +δτ/2 through Krylov exponentiation
# of the local effective map, factorized by QR (LQ), and the bond matrix C is
# evolved by -δτ/2; the last (first) site performs a full +δτ step.
#
# The initial identity state is zero-padded up to the bond dimension D and
# canonicalized without truncation, so that the zero-weight directions become
# orthonormal directions of the environments and the sweeps can populate the
# full bond profile min(d^j, d^{L-j}, D) as the flow builds up correlations.
# Truncation is applied only in the final canonicalization.


# on real-time lattices `influenceoperator` returns 4 branch MPOs
# ((+,+), (+,−), (−,+), (−,−)); the total influence operator driving the
# flow is their SUM (MPS/MPO direct sum, bond dimensions add up). Indeed the
# differential IF built by `differentialinfluencefunctional` is the Hadamard
# (site-wise) product e^{dt·h1}∘e^{dt·h2}∘e^{dt·h3}∘e^{dt·h4}, and in the
# element-wise algebra of ADT/PT products e^a∘e^b = e^{a+b}, so the generator
# of the full IF is h1+h2+h3+h4 (NOT the Hadamard product h1∘h2∘h3∘h4).
# NOTE: the direct-sum bond dimensions of the four branches add up, so the
# branches are summed one by one, compressing with SVD canonicalization after
# each addition to keep the bond dimension bounded. The compression uses the
# tight `DefaultMPOTruncation` (not `alg.trunc`): the compression error of H
# is exponentially amplified by the flow (IF = e^H), so a loose tolerance
# would degrade the accuracy of the influence functional.
function _tdvpif_hamiltonian(lattice, corr, hyb, alg::TDVPIF)
	h1, h2, h3, h4 = influenceoperator(lattice, corr, hyb, algexpan=alg.algexpan)
	orth = Orthogonalize(SVD(), DefaultMPOTruncation; normalize=false)
	H = h1 + h2
	canonicalize!(H, alg=orth)
	H = H + h3
	canonicalize!(H, alg=orth)
	H = H + h4
	canonicalize!(H, alg=orth)
	return H
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
# preparation: lift z to the flow bond dimension (zero-padded) and canonicalize
# without truncation, so that the zero-weight directions become orthonormal
# directions of the environments and the sweeps can populate the full bond
# profile min(d^j, d^{L-j}, D); finalization: canonicalize with the truncation
# scheme. The global scaling factor of z is carried through the flow by the
# `_renormalize!` bookkeeping, so the output value is e^H·z(0) regardless of
# the gauge of the input.
function _tdvpif_hybriddynamics_adt!(z::ADT, H::ADT, alg::TDVPIF)
	increase_bond!(z, alg.trunc.D)
	canonicalize!(z, alg=Orthogonalize(SVD(), NoTruncation(); normalize=false))
	_tdvpif_flow_adt!(z, H, alg)
	canonicalize!(z, alg=Orthogonalize(SVD(), alg.trunc; normalize=false))
	alg.callback(Float64[])
	return z
end

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

# effective map on the center tensor at site j: the tangent-space projection
# of the product H·z, built from the left environment hleft::(bra, H, ket),
# the influence operator H[j] and the right environment hright::(ket, H, bra)
function _tdvpif_ac_prime_adt(AC::DenseMPSTensor, Hj::DenseMPSTensor, hleft::DenseMPSTensor, hright::DenseMPSTensor)
	left_xy = get_left_xy(hleft, Hj, AC)
	@tensor tmp[1, 2, 5] := left_xy[1, 2, 3, 4] * hright[4, 3, 5]
	return tmp
end

# effective map on the bond matrix C: contraction of the left and right
# environments with the H bonds passing straight through the bond
function _tdvpif_c_prime_adt(C::AbstractMatrix, hleft::DenseMPSTensor, hright::DenseMPSTensor)
	@tensor tmp[1, 6] := (hleft[1, 2, 3] * C[3, 4]) * hright[4, 2, 6]
	return tmp
end

# initialize the right environments ⟨z|H|z⟩; hstorage[i]::(ket, H, bra) is the
# partial contraction over sites i:L
function _tdvpif_init_hstorage_adt(z::ADT, H::ADT)
	L = length(z)
	T = scalartype(z)
	hstorage = Vector{Array{T, 3}}(undef, L+1)
	hstorage[L+1] = ones(T, space_r(z), space_r(H), space_r(z))
	for i in L:-1:2
		hstorage[i] = updatemultright(hstorage[i+1], z[i], H[i], z[i])
	end
	hstorage[1] = ones(T, space_l(z), space_l(H), space_l(z))
	return hstorage
end

function _tdvpif_leftsweep_adt!(z::ADT, H::ADT, hstorage, δτ::Float64, alg::TDVPIF)
	L = length(z)
	krylov = Arnoldi()
	for site in 1:L-1
		(alg.verbosity > 3) && println("TDVPIF: left sweep at site $site")
		# forward half-step of the center tensor
		AC, info = exponentiate(x -> _tdvpif_ac_prime_adt(x, H[site], hstorage[site], hstorage[site+1]), δτ/2, z[site], krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at site $site"
		AL, C = leftorth!(AC, (1, 2), (3,))
		z[site] = AL
		C = Matrix(C)
		_renormalize!(z, C, false)
		# left environment with the new AL on both the bra and ket sides
		hnew = updatemultleft(hstorage[site], AL, H[site], AL)
		# backward half-step of the bond matrix
		C, info = exponentiate(x -> _tdvpif_c_prime_adt(x, hnew, hstorage[site+1]), -δτ/2, C, krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at bond $(site+1)"
		_renormalize!(z, C, false)
		hstorage[site+1] = hnew
		# absorb the bond matrix into the next site
		z[site+1] = @tensor tmp[-1, -2, -3] := C[-1, 1] * z[site+1][1, -2, -3]
	end
	# full step at the last site
	AC, info = exponentiate(x -> _tdvpif_ac_prime_adt(x, H[L], hstorage[L], hstorage[L+1]), δτ, z[L], krylov)
	(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at site $L"
	z[L] = AC
	_renormalize!(z, z[L], false)
	return z
end

function _tdvpif_rightsweep_adt!(z::ADT, H::ADT, hstorage, δτ::Float64, alg::TDVPIF)
	krylov = Arnoldi()
	for site in length(z)-1:-1:1
		(alg.verbosity > 3) && println("TDVPIF: right sweep at site $site")
		C, AR = rightorth!(z[site+1], (1,), (2, 3))
		z[site+1] = AR
		C = Matrix(C)
		# right environment with the new AR on both the bra and ket sides
		hnew = updatemultright(hstorage[site+2], AR, H[site+1], AR)
		# backward half-step of the bond matrix
		C, info = exponentiate(x -> _tdvpif_c_prime_adt(x, hstorage[site+1], hnew), -δτ/2, C, krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at bond $(site+1)"
		_renormalize!(z, C, false)
		hstorage[site+1] = hnew
		# absorb the bond matrix into the site on the left
		z[site] = @tensor tmp[-1, -2, -3] := z[site][-1, -2, 1] * C[1, -3]
		# forward half-step of the center tensor
		AC, info = exponentiate(x -> _tdvpif_ac_prime_adt(x, H[site], hstorage[site], hstorage[site+1]), δτ/2, z[site], krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at site $site"
		z[site] = AC
		_renormalize!(z, z[site], false)
	end
	return z
end

function _tdvpif_flow_adt!(z::ADT, H::ADT, alg::TDVPIF)
	# The flow contracts the raw site tensors of `H` and never applies its
	# global scaling factor. The ADT convention is
	#     value = scaling^L · ∏(site tensors),
	# so any H with scaling ≠ 1 (e.g. after a canonicalization, which
	# redistributes local weights into the scaling factor) would silently be
	# evolved as H / scaling^L. Absorb the scaling into the site tensors first
	# to make the flow independent of the gauge in which H is represented.
	_absorb_scaling!(H)
	hstorage = _tdvpif_init_hstorage_adt(z, H)
	nsteps = round(Int, 1 / alg.δ)
	δτ = 1 / nsteps
	for n in 1:nsteps
		_tdvpif_leftsweep_adt!(z, H, hstorage, δτ, alg)
		_tdvpif_rightsweep_adt!(z, H, hstorage, δτ, alg)
		(alg.verbosity > 1) && println("TDVPIF step $n/$nsteps, τ = $(n * δτ)")
	end
	(alg.verbosity > 1) && println("TDVPIF flow finished, τ = 1")
	return z
end

"""
	hybriddynamics(lattice::ImagADTLattice1Order, corr::ImagCorrelationFunction, hyb::AdditiveHyb, alg::TDVPIF)

Construct the influence functional on an imaginary-time ADT lattice with the `TDVPIF` algorithm: evolve the identity influence functional along the imaginary-time flow dz/dτ = H·z (H = influence operator) from τ = 0 to τ = 1 with second-order single-site TDVP sweeps.

# Returns
The influence functional, represented as an `ADT`.

See also [`TDVPIF`](@ref).
"""
function hybriddynamics(lattice::ImagADTLattice1Order, corr::ImagCorrelationFunction, hyb::AdditiveHyb, alg::TDVPIF)
	T = promote_type(scalartype(corr), scalartype(hyb), Float64)
	z = vacuumstate(T, lattice)
	return hybriddynamics!(z, lattice, corr, hyb, alg)
end

"""
	hybriddynamics!(gmps::ADT, lattice::ImagADTLattice1Order, corr::ImagCorrelationFunction, hyb::AdditiveHyb, alg::TDVPIF)

In-place version of the `TDVPIF` algorithm on ADT lattices: the influence operator H drives the TDVP flow directly on `gmps`, i.e. the flow evolves `z(τ=1) = e^H·gmps` from `z(0) = gmps` in place. This allows merging the influence functional into an arbitrary initial MPO — e.g. the pure impurity dynamics (with Lindblad dissipation) obtained from `sysdynamics` — in a single flow instead of constructing the IF separately and multiplying it in afterwards.

# Returns
The modified `gmps`.
"""
function hybriddynamics!(gmps::ADT, lattice::ImagADTLattice1Order, corr::ImagCorrelationFunction, hyb::AdditiveHyb, alg::TDVPIF)
	# influence operator H in MPO (fused-leg MPS) form
	H = influenceoperator(lattice, corr, hyb, algexpan=alg.algexpan)
	return _tdvpif_hybriddynamics_adt!(gmps, H, alg)
end

"""
	hybriddynamics(lattice::RealADTLattice1Order, corr::RealCorrelationFunction, hyb::AdditiveHyb, alg::TDVPIF)

Construct the influence functional on a real-time ADT lattice with the `TDVPIF` algorithm: the 4 branch influence operators returned by `influenceoperator` ((+,+), (+,−), (−,+), (−,−)) are summed one by one into a single influence operator H (each partial sum compressed by SVD canonicalization with `DefaultMPOTruncation`), which then drives the same TDVP imaginary-time flow as in the imaginary-time case.

# Returns
The influence functional, represented as an `ADT`.

See also [`TDVPIF`](@ref).
"""
function hybriddynamics(lattice::RealADTLattice1Order, corr::RealCorrelationFunction, hyb::AdditiveHyb, alg::TDVPIF)
	T = promote_type(scalartype(corr), scalartype(hyb), Float64)
	z = vacuumstate(T, lattice)
	return hybriddynamics!(z, lattice, corr, hyb, alg)
end

"""
	hybriddynamics!(gmps::ADT, lattice::RealADTLattice1Order, corr::RealCorrelationFunction, hyb::AdditiveHyb, alg::TDVPIF)

In-place version of the `TDVPIF` algorithm on real-time ADT lattices: the influence operator H (the sum of the 4 branch operators, see [`hybriddynamics`](@ref)) drives the TDVP flow directly on `gmps`, i.e. the flow evolves `z(τ=1) = e^H·gmps` from `z(0) = gmps` in place, merging the influence functional into an arbitrary initial MPO in a single flow.

# Returns
The modified `gmps`.
"""
function hybriddynamics!(gmps::ADT, lattice::RealADTLattice1Order, corr::RealCorrelationFunction, hyb::AdditiveHyb, alg::TDVPIF)
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
function _tdvpif_hybriddynamics_pt!(z::ProcessTensor, H::ProcessTensor, alg::TDVPIF)
	increase_bond!(z, alg.trunc.D)
	canonicalize!(z, alg=Orthogonalize(SVD(), NoTruncation(); normalize=false))
	_tdvpif_flow_pt!(z, H, alg)
	canonicalize!(z, alg=Orthogonalize(SVD(), alg.trunc; normalize=false))
	alg.callback(Float64[])
	return z
end

# effective map on the center tensor: the tangent-space projection of the
# operator product H·z; cleft/cright::(bra, H, ket)
function _tdvpif_ac_prime_pt(AC::DenseMPOTensor, Hj::DenseMPOTensor, cleft::DenseMPSTensor, cright::DenseMPSTensor)
	return reduceH_single_site(AC, Hj, cleft, cright)
end

# effective map on the bond matrix C: contraction of the left and right
# environments with the H bonds passing straight through the bond
function _tdvpif_c_prime_pt(C::AbstractMatrix, cleft::DenseMPSTensor, cright::DenseMPSTensor)
	@tensor tmp[1, 6] := (cleft[1, 2, 3] * C[3, 4]) * cright[6, 2, 4]
	return tmp
end

function _tdvpif_leftsweep_pt!(z::ProcessTensor, H::ProcessTensor, hstorage, δτ::Float64, alg::TDVPIF)
	L = length(z)
	krylov = Arnoldi()
	for site in 1:L-1
		(alg.verbosity > 3) && println("TDVPIF: left sweep at site $site")
		# forward half-step of the center tensor
		AC, info = exponentiate(x -> _tdvpif_ac_prime_pt(x, H[site], hstorage[site], hstorage[site+1]), δτ/2, z[site], krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at site $site"
		AL, C = leftorth!(AC, (1, 2, 4), (3,))
		AL = permute(AL, (1, 2, 4, 3))
		z[site] = AL
		C = Matrix(C)
		_renormalize!(z, C, false)
		# left environment with the new AL on both the bra and ket sides
		hnew = updateleft(hstorage[site], AL, H[site], AL)
		# backward half-step of the bond matrix
		C, info = exponentiate(x -> _tdvpif_c_prime_pt(x, hnew, hstorage[site+1]), -δτ/2, C, krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at bond $(site+1)"
		_renormalize!(z, C, false)
		hstorage[site+1] = hnew
		# absorb the bond matrix into the next site
		z[site+1] = @tensor tmp[-1, -2, -3, -4] := C[-1, 1] * z[site+1][1, -2, -3, -4]
	end
	# full step at the last site
	AC, info = exponentiate(x -> _tdvpif_ac_prime_pt(x, H[L], hstorage[L], hstorage[L+1]), δτ, z[L], krylov)
	(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at site $L"
	z[L] = AC
	_renormalize!(z, z[L], false)
	return z
end

function _tdvpif_rightsweep_pt!(z::ProcessTensor, H::ProcessTensor, hstorage, δτ::Float64, alg::TDVPIF)
	krylov = Arnoldi()
	for site in length(z)-1:-1:1
		(alg.verbosity > 3) && println("TDVPIF: right sweep at site $site")
		C, AR = rightorth!(z[site+1], (1,), (2, 3, 4))
		z[site+1] = AR
		C = Matrix(C)
		# right environment with the new AR on both the bra and ket sides
		hnew = updateright(hstorage[site+2], AR, H[site+1], AR)
		# backward half-step of the bond matrix
		C, info = exponentiate(x -> _tdvpif_c_prime_pt(x, hstorage[site+1], hnew), -δτ/2, C, krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at bond $(site+1)"
		_renormalize!(z, C, false)
		hstorage[site+1] = hnew
		# absorb the bond matrix into the site on the left
		z[site] = @tensor tmp[-1, -2, -3, -4] := z[site][-1, -2, 1, -4] * C[1, -3]
		# forward half-step of the center tensor
		AC, info = exponentiate(x -> _tdvpif_ac_prime_pt(x, H[site], hstorage[site], hstorage[site+1]), δτ/2, z[site], krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at site $site"
		z[site] = AC
		_renormalize!(z, z[site], false)
	end
	return z
end

function _tdvpif_flow_pt!(z::ProcessTensor, H::ProcessTensor, alg::TDVPIF)
	# same scaling absorption as in `_tdvpif_flow_adt!`
	_absorb_scaling!(H)
	hstorage = init_hstorage_right(z, H, z)
	nsteps = round(Int, 1 / alg.δ)
	δτ = 1 / nsteps
	for n in 1:nsteps
		_tdvpif_leftsweep_pt!(z, H, hstorage, δτ, alg)
		_tdvpif_rightsweep_pt!(z, H, hstorage, δτ, alg)
		(alg.verbosity > 1) && println("TDVPIF step $n/$nsteps, τ = $(n * δτ)")
	end
	(alg.verbosity > 1) && println("TDVPIF flow finished, τ = 1")
	return z
end

"""
	hybriddynamics(lattice::ImagPTLattice1Order, corr::ImagCorrelationFunction, hyb::GeneralHybStyle, alg::TDVPIF)

Construct the influence functional on an imaginary-time PT lattice with the `TDVPIF` algorithm, supporting `GeneralHybStyle` (e.g. `NonAdditiveHyb`, `NonDiagonalHyb`) coupling.

# Returns
The influence functional, represented as a `ProcessTensor`.

See also [`TDVPIF`](@ref).
"""
function hybriddynamics(lattice::ImagPTLattice1Order, corr::ImagCorrelationFunction, hyb::GeneralHybStyle, alg::TDVPIF)
	T = promote_type(scalartype(corr), scalartype(hyb), Float64)
	z = vacuumstate(T, lattice)
	return hybriddynamics!(z, lattice, corr, hyb, alg)
end

"""
	hybriddynamics!(gmps::ProcessTensor, lattice::ImagPTLattice1Order, corr::ImagCorrelationFunction, hyb::GeneralHybStyle, alg::TDVPIF)

In-place version of the `TDVPIF` algorithm on PT lattices: the influence operator H drives the TDVP flow directly on `gmps`, i.e. the flow evolves `z(τ=1) = e^H·gmps` from `z(0) = gmps` in place. This allows merging the influence functional into an arbitrary initial MPO — e.g. the pure impurity dynamics (with Lindblad dissipation) obtained from `sysdynamics` — in a single flow instead of constructing the IF separately and multiplying it in afterwards.

# Returns
The modified `gmps`.
"""
function hybriddynamics!(gmps::ProcessTensor, lattice::ImagPTLattice1Order, corr::ImagCorrelationFunction, hyb::GeneralHybStyle, alg::TDVPIF)
	# influence operator H in MPO form
	H = influenceoperator(lattice, corr, hyb, algexpan=alg.algexpan)
	return _tdvpif_hybriddynamics_pt!(gmps, H, alg)
end

"""
	hybriddynamics(lattice::RealPTLattice1Order, corr::RealCorrelationFunction, hyb::GeneralHybStyle, alg::TDVPIF)

Construct the influence functional on a real-time PT lattice with the `TDVPIF` algorithm: the 4 branch influence operators returned by `influenceoperator` ((+,+), (+,−), (−,+), (−,−)) are summed one by one into a single influence operator H (each partial sum compressed by SVD canonicalization with `DefaultMPOTruncation`), which then drives the same TDVP imaginary-time flow as in the imaginary-time case.

# Returns
The influence functional, represented as a `ProcessTensor`.

See also [`TDVPIF`](@ref).
"""
function hybriddynamics(lattice::RealPTLattice1Order, corr::RealCorrelationFunction, hyb::GeneralHybStyle, alg::TDVPIF)
	T = promote_type(scalartype(corr), scalartype(hyb), Float64)
	z = vacuumstate(T, lattice)
	return hybriddynamics!(z, lattice, corr, hyb, alg)
end

"""
	hybriddynamics!(gmps::ProcessTensor, lattice::RealPTLattice1Order, corr::RealCorrelationFunction, hyb::GeneralHybStyle, alg::TDVPIF)

In-place version of the `TDVPIF` algorithm on real-time PT lattices: the influence operator H (the sum of the 4 branch operators, see [`hybriddynamics`](@ref)) drives the TDVP flow directly on `gmps`, i.e. the flow evolves `z(τ=1) = e^H·gmps` from `z(0) = gmps` in place, merging the influence functional into an arbitrary initial MPO in a single flow.

# Returns
The modified `gmps`.
"""
function hybriddynamics!(gmps::ProcessTensor, lattice::RealPTLattice1Order, corr::RealCorrelationFunction, hyb::GeneralHybStyle, alg::TDVPIF)
	H = _tdvpif_hamiltonian(lattice, corr, hyb, alg)
	return _tdvpif_hybriddynamics_pt!(gmps, H, alg)
end
