# Practice guide

This page summarizes hands-on experience accumulated while reproducing published examples (spinboson, strathearn2018, otterpohl2025, and guo2026 under `docs/tutorials/`): how to choose a computational path, how to run convergence checks, common pitfalls, and how to diagnose wrong results. **We recommend reading it through once before running any production calculation.**

## 1. Choosing a computational path

Answer two questions first:

1. **What is the form of the coupling operator?** Diagonal (Hermitian, commuting with $H_S$) → the ADT framework; off-diagonal (e.g. $\sigma_-$) or multiple baths with non-commuting couplings → the PT framework.
2. **Which observables do you need?** Diagonal ones only (e.g. $\langle S_z(t)\rangle$) → ADT suffices; off-diagonal ones (e.g. $\langle\sigma_x(t)\rangle$, $\langle a\rangle$) → you must take the PT path.

| Situation | Recommended path | Rationale |
|---|---|---|
| Diagonal coupling + diagonal observables + d=2 | `ADTLattice` + `AdditiveHyb` + `PartialIF` (default) | Fastest; bond dimension 2 per factor |
| Diagonal coupling + diagonal observables + d≥3 (e.g. spin-1) | `ADTLattice` + `AdditiveHyb` + `XTRGIF` | The number of PartialIF factors grows with d; XTRG-IF is clearly faster |
| Off-diagonal observables needed | `PTLattice` + `NonDiagonalHyb` + `XTRGIF`; or the ADT path + operator insertion (see §3) | The PT path is more efficient for bulk measurements; ADT operator insertion suits a few operators / two-point correlations |
| Off-diagonal / non-Hermitian coupling | `PTLattice` + `NonDiagonalHyb` + `XTRGIF` | The only supported option |

!!! tip "Off-diagonal observables are possible even with diagonal coupling"
    Two routes: ① for a **real** diagonal operator `op`, `NonDiagonalHyb(op)` is physically equivalent to `AdditiveHyb(op)`, but it takes the PT path and allows measuring off-diagonal observables in bulk with `ContourOperator` (e.g. in the otterpohl2025 tutorial, `NonDiagonalHyb(Matrix(Diagonal(0.5 .* [1.0, -1.0])))` is the $\sigma_z/2$ coupling, measuring both $P=\langle\sigma_z\rangle$ and $C=\langle\sigma_x\rangle$ at once); ② stay on the ADT path and use **operator insertion** (`sysdynamics(lattice, model, ct)`, see §3), which suits a few operators or two-point correlation functions.

## 2. Convergence checks: workflow and rules of thumb

Four sources of error: `δt`, `χ`, `n` (number of exponential-expansion terms), and `k` (number of XTRG steps). Recommended workflow:

1. Fix a large `χ`, `n=20`, `k=7`, and scan `δt` (e.g. 0.2 → 0.1 → 0.05); take a `δt` beyond which the results no longer change;
2. Fix that `δt` and scan `χ` (e.g. 25 → 50 → 100) until the observables (especially peak and minimum positions) no longer change;
3. If XTRG-IF is used, check the expansion error of `exponential_expansion` (see §5).

Rules of thumb from practice (use them as starting points; still verify with your own parameters):

| Example | Difficulty | Rules of thumb |
|---|---|---|
| Ohmic SBM ($\Omega=1$, $\omega_c=5$, $T=0$) **crossing the localization transition** | Entanglement is largest near the transition point | $\chi\ge 50$ and $\delta t\le 0.1$ (at $\chi=25$, $\delta t=0.15$ the $S_z$ error can reach ~0.04) |
| Low-temperature Ohmic SBM | — | $\beta=20$ is effectively zero temperature (differs from $\beta=100$ by <0.001) |
| $1/f^\eta$ noise ($s<0$, Otterpohl-type) | Power-law memory in time keeps growing the IF bond dimension | For short times ($t\le5$), $\chi=40$ suffices; the stronger coupling $\alpha=0.04$ needs $\chi=100$ |
| JC-type (number-conserving) off-diagonal coupling | Far less entanglement than the Rabi type | $\chi=30$ usually suffices |
| spin-1 (d=3) + XTRG-IF | Large local dimension | $\chi=50$, $\delta t=0.1$, $n=20$, $k=7$ |

## 3. Constructing observables: the most common pitfalls

### Operator arguments of the hybridization styles

- `AdditiveHyb` accepts a **vector of diagonal elements** or a diagonal matrix. Passing a matrix `z` works (its diagonal is taken), but for clarity pass `zdiag = diag(z)` directly;
- `NonDiagonalHyb(op)` means $\hat{A}\hat{a} + \hat{A}^\dagger\hat{a}^\dagger$, i.e. `op` is the coupling operator $\hat{A}$ itself (**including its coefficient**). The spin-boson coupling $\tfrac{1}{2}\sigma_z\hat\xi$ should be written `NonDiagonalHyb(Matrix(Diagonal(0.5 .* [1.0, -1.0])))`, i.e. $\hat{A}=\sigma_z/2$;
- Beware of operator convention differences across the literature: some papers use $S=\sigma/2$ (the TEMPO package convention, with no factor of two on the Pauli operators). When comparing against published data, first check the definitions of the coupling operator and $H_S$ (e.g. whether $H_s=\Omega\sigma_x/2$ corresponds to `ImpurityHamiltonian(Ω .* Sx/2)` or `Ω .* Sx` depends on whether $S$ or $\sigma$ is used).

### ADT path: `ADTTerm` is a single-point diagonal operator; use operator insertion for non-diagonal quantities

```julia
m = ADTTerm(index(lattice, i, branch=:+), zdiag)   # zdiag is the vector of diagonal elements
v = expectationvalue(m, cache)                     # already normalized
```

- In `ADTTerm(pos, data)`, `data` is a **vector of diagonal elements** acting on a single site (the system site and its branch-paired site are absorbed by the lattice into a single step);
- The multi-site form `ADTTerm((pos2, pos1), (v2, v1))` can measure **diagonal** two-point correlations: `mps2 = apply!(m, copy(mps)); v = integrate(mps2) / integrate(mps)`;
- **Do not** pass a full matrix to a single-point `ADTTerm`;
- For **off-diagonal single-point operators or arbitrary two-point correlations**, use operator insertion (works for both ADT and PT; see `test/models/rabimodel.jl`):

  ```julia
  ct   = ContourOperator([ContourIndex(i), ContourIndex(1)], [op2, op1])  # op can be an arbitrary matrix
  mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
  mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)   # real time needs ρ₀; imaginary time does not
  mps2 = mult!(mpsK, mpsI, trunc=trunc)
  v    = integrate(mps2) / integrate(mps)             # already normalized
  ```

  Each insertion rebuilds `mpsK` and re-multiplies the IF, so it suits a few operators; for bulk single-point measurements prefer the environment-cache path.

### PT path: `ContourOperator` accepts arbitrary matrices

```julia
cache = environments(lattice, mps, ρ₀=ρimp)        # real-time PT requires ρ₀!
m = ContourOperator(ContourIndex(i, :+), x)        # x is an arbitrary matrix
v = expectationvalue(m, cache)
```

- In real-time PT, `environments` requires `ρ₀` (the initial density matrix); omitting it silently yields results for the default maximally mixed state;
- Two-point insertions use `ContourOperator([c1, c2], [op1, op2])` and enter the system dynamics via `sysdynamics(lattice, model, ct)` (the standard approach for imaginary-time Green functions).

### Initial states

```julia
mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)  # ADT path
```

- `ρ₀` can be a density matrix or a diagonal vector (the diagonal elements of a pure state);
- The initial state is applied to `mpsK` (the system dynamics), not to `mpsI` (the influence functional).

## 4. Unphysical results? Suspect insufficient χ first

The **typical symptom of insufficient χ is not an error message but silently wrong results**:

- Normalization drift: $P(t)>1$, $|C(t)|>1$;
- Unphysical oscillations or clearly shifted minimum positions;
- Diverging long-time behavior.

For example, in the otterpohl2025 tutorial, $\alpha=0.04$ with $\chi=40$ gives $P(\text{end})=2.94$ and $C_{\min}=-1.58$ (unphysical); switching to $\chi=100$ restores $P\le1$ and physical oscillations in $C$. **Every new example should be checked with at least two χ values.**

## 5. Spectral densities and correlation functions

- `spectrum(f, lb=0, ub=...)`: the upper bound should be large enough that the spectral weight has essentially decayed (for Ohmic/sub-Ohmic spectra, `ub=3~5ωc` is enough);
- **Low-frequency divergent spectra**: if $J(0)\neq 0$ and the spectrum is convolved with the Bose factor $1/\epsilon$, it diverges logarithmically at $\omega=0$. The remedy is to integrate from a small nonzero lower bound (e.g. the segal2023 tutorial starts at $\gamma/100$);
- **Zero temperature**: `bosonicbath(spec, β=Inf)`;
- **`ExpExp` warning** `can not find a good approximation with L=..., n=...`: the Prony expansion did not reach `tol`. This is common for power-law memory kernels with $s<0$; the computation is usually still fine (TEMPO only needs the correlation function), but you should:
  - increase `n` (e.g. 20 → 30) or relax `tol` (e.g. `1e-8` → `1e-6`) and compare the results;
  - treat the warning seriously only during convergence checks; final results should be re-verified with parameters that produce no warning (or a negligible error).

## 6. Performance: caching and algorithm selection

- **Cache the influence-functional MPS/MPO**: the XTRG construction of XTRG-IF can take hours (e.g. about 3.7 hours for $\chi=100$, $N=103$). Serialize `mpsI` to disk with `Serialization` and deserialize it directly when re-running observables:

  ```julia
  using Serialization
  mpspath = "data/mypi_beta20_dt0.1_alpha0.3_chi50_N200.mps"
  mpsI = ispath(mpspath) ? Serialization.deserialize(mpspath) :
         (I = hybriddynamics(lattice, corr, hyb, alg); Serialization.serialize(mpspath, I); I)
  ```
- **Two switches for XTRG-IF**: `fast=true` (tree bisection, about `k` multiplications instead of $2^k-1$) and a sensible `k` (`k=7` means a smallest step of $1/128$; larger `k` has diminishing returns);
- **PartialIF is the fastest for d=2**; do not switch to XTRG-IF blindly just because it is the "new algorithm". For d≥3, XTRG-IF has a clear advantage;
- Set `OMP_NUM_THREADS=1` (or a suitable number of threads for your machine) to avoid BLAS oversubscription; the DMRG multiplication itself is single-threaded.

## 7. Comparison with reference data (papers / ED)

1. **Check conventions first**: the coupling operator, $H_S$, the normalization of the spectral density (e.g. whether $J(\omega)=2\alpha\omega e^{-\omega/\omega_c}$ carries a factor of $\pi$), and the time units;
2. **Align grids by time**: TEMPO's sample points are `t_i = i*δt, i=1..N`; when comparing against external data (ED grids, digitized data), index by time value rather than by array position:
   ```julia
   n_on_t = [ne[round(Int, t_i / 0.005) + 1] for t_i in t]   # ED step size 0.005
   ```
3. **Digitized paper data**: check that the in-figure legend is consistent with the paper's text (in the strathearn2018 tutorial, α was once mislabeled because the legend of another panel was used by mistake); digitized curves carry noise, so focus the comparison on trends and extremum positions rather than pointwise differences;
4. **Physical dimensions**: the bath UV cutoff $1/\omega_c$ sets characteristic time scales such as the position of the minimum in the pseudocoherent regime; keep the units self-consistent.

## 8. Quick diagnostics table

| Symptom | Likely cause | Fix |
|---|---|---|
| $P(t)>1$, $\lvert C\rvert>1$ | Insufficient χ | Increase χ and rerun (§4) |
| Results diverge with t | Insufficient χ or too-large δt | First double χ, then halve δt |
| Systematic shift of minimum/peak positions | δt discretization error | Scan δt (§2) |
| `ExpExp` non-convergence warning | Correlation function hard to expand | Increase n or relax tol, and re-verify (§5) |
| Expectation values correspond to the maximally mixed state | Forgot to pass `ρ₀` on the PT path | `environments(lattice, mps, ρ₀=ρimp)` (§3) |
| `expectationvalue` type error | ADT/PT terms mixed up | Pair `ADTTerm` with the ADT cache and `ContourOperator` with the PT cache (§3) |
| Bath correlation function diverges | Low-frequency divergence of the spectrum | Use a small positive lower bound for the spectral integral (§5) |
| Disagreement with the paper | Operator/spectral convention differences | Check $S$ vs $\sigma$, $\pi$ normalization, and units (§7) |
