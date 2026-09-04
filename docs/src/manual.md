# Manual

This page introduces the usage and key points of the core API, organized by component; full per-symbol documentation can be found in the [API reference](@ref). Contents:

1. [Truncation schemes `TruncationScheme`](@ref)
2. [Lattices `ADTLattice` / `PTLattice`](@ref)
3. [System operators `AbstractImpurityOperator`](@ref)
4. [Hybridization styles `HybridizationStyle`](@ref)
5. [Influence functional algorithms `InfluenceFunctionalAlgorithm`](@ref)
6. [Baths and correlation functions (`ImpurityModelBase`)](@ref)
7. [ADT / PT tensor operations](@ref)
8. [Observables](@ref manual_observables)
9. [MPO Hamiltonians (long-range interaction tools)](@ref)
10. [Hyperparameters and sources of error](@ref)
11. [Code structure](@ref)
12. [Correspondence with the literature](@ref)

## Truncation schemes `TruncationScheme`

Truncation schemes for compressing the bond dimension in tensor network computations:

| Type | Construction | Description |
|---|---|---|
| `TruncationDimCutoff` | `truncdimcutoff(D=χ, ϵ=ε, add_back=0)` | Limits both the maximum bond dimension `D` and the truncation threshold `ϵ` (recommended) |
| `TruncateDim` | `truncdim(D)` | Limits only the bond dimension |
| `TruncateCutoff` | `trunccutoff(ϵ=ε)` | Truncates only by the singular-value threshold |
| `NoTruncation` | `NoTruncation()` | No truncation |

Predefined defaults (`src/defaults.jl`):

```julia
DefaultTruncation        # D=100, ϵ=1e-14
DefaultITruncation       # D=200,  ϵ=1e-10   (default for building the IF)
DefaultKTruncation       # D=1000, ϵ=1e-10   (default for system dynamics)
DefaultIntegrationTruncation  # D=10000, ϵ=1e-12
DefaultMPOTruncation     # D=10000, ϵ=1e-12
```

## Lattices `ADTLattice` / `PTLattice`

Unified construction entry points:

```julia
lattice = ADTLattice(N=Nt, δt=δt, contour=:real)      # real time
lattice = PTLattice(N=N,  δτ=δτ, contour=:imag)       # imaginary time
lattice = PTLattice(Nt=Nt, Nτ=Nτ, δt=δt, δτ=δτ, d=d, contour=:mixed)
```

Commonly used properties and functions:

- `length(lattice)`: total number of sites; `phydim(lattice)`: local dimension `d`;
- `lattice.d`, `lattice.N` / `Nt` / `Nτ`, `lattice.δt` / `δτ`, `lattice.t`, `lattice.β` (`N*δτ` for imaginary time);
- `branches(lattice)`: returns the branches, e.g. `(:+, :-)` (real time) or `(:τ,)` (imaginary time);
- `index(lattice, i, branch=:+)`: maps a (time step, branch) pair to a site position;
- `ContourIndex(i, :+)`: a contour index object;
- `vacuumstate(lattice)`: initializes the vacuum-state ADT/PT.

!!! note "Mind the grid endpoints"
    `N` time steps cover the time interval `t = N*δt`, with sample points `i = 1..N` corresponding to `t_i = i*δt`. When comparing against external data (e.g. exact diagonalization), check whether the grid includes the endpoints and the starting point; it is recommended to align explicitly in time rather than by index.

## System operators `AbstractImpurityOperator`

| Type | Construction | Description |
|---|---|---|
| `ImpurityHamiltonian` | `ImpurityHamiltonian(m::Matrix)` | Unitary evolution $e^{\mp i \hat{H}_S \delta t}$ |
| `ImpurityLindbladian` | `ImpurityLindbladian(H, jumpops)` | Dissipative (Lindblad) evolution; `jumpops` is a list of jump operators |

`sysdynamics(lattice, model, trunc=...)` generates the bare system dynamics tensor `mpsK` (ADT or PT).

## Hybridization styles `HybridizationStyle`

Describes the specific form of the system–bath coupling:

| Type | Coupling form | Constraint | Framework |
|---|---|---|---|
| `AdditiveHyb(op)` | $\hat{A}(\hat{b}^\dagger+\hat{b})$ type, diagonal | `op` is a vector or a diagonal matrix (Hermitian) | ADT |
| `NonAdditiveHyb(op)` | $\hat{A}(\hat{a}+\hat{a}^\dagger)$, $\hat{A}=\hat{A}^\dagger$ | `op` is a Hermitian matrix | PT |
| `NonDiagonalHyb(op)` | $\hat{A}\hat{a} + \hat{A}^\dagger \hat{a}^\dagger$ | `op` is any square matrix | PT (the core case in the paper) |

Time-dependent versions (`TdHybridizationStyle`): `AdditiveTdHyb`, `NonAdditiveTdHyb`, `NonDiagonalTdHyb`.

## Influence functional algorithms `InfluenceFunctionalAlgorithm`

```julia
# 1) Partial influence functional (Partial IF): only for the diagonal coupling AdditiveHyb
alg = PartialIF(trunc=trunc)

# 2) Translationally invariant influence functional (XTRG-IF, XTRG-style), recommended by the paper
alg = XTRGIF(;
    algexpan = OverDeterminedProny(n=20, tol=1.0e-8, verbosity=2),  # exponential expansion of the hybridization function
    algevo   = WII(),        # or WI(), ComplexStepper(), FirstOrderStepper()
    algmult  = DMRGMult1(trunc, initguess=:rand, maxiter=10),   # or SVDCompression(trunc)
    k        = 7,            # number of XTRG steps: time step 1/2^k
    fast     = true,         # true: tree/bisection scheme (about k multiplications); false: serial 2^k-1 multiplications
    verbosity= 0,
)
```

- `algexpan`: an `ExponentialExpansionAlgorithm`, including `OverDeterminedProny` and `DeterminedProny` (`exponential_expansion` and `expansion_error` can be used for error analysis; provided by the re-exported `ExpExp` package);
- `algevo`: a `TimeEvoMPOAlgorithm`, the stepper that exponentiates $\hat{H}_{\text{eff}}$ into an MPO (`WI`/`WII`/`FirstOrderStepper`/`ComplexStepper`);
- `algmult`: a `DMRGAlgorithm` that compresses MPO–MPO multiplications, either `DMRGMult1` (single-site DMRG iteration, `initguess ∈ {:svd, :pre, :rand}`, `maxiter`) or `SVDCompression`.

Entry-point functions for building the IF:

```julia
mpsI = hybriddynamics(lattice, corr, hyb)                 # default PartialIF / default truncation
mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)    # diagonal coupling + partial IF
mpsI = hybriddynamics(lattice, corr, hyb, alg)            # specify the IF algorithm (diagonal or off-diagonal)
mpsI = hybriddynamics!(gmps, lattice, corr, hyb, alg)     # in-place version
mpsI = hybriddynamics_naive(lattice, corr, hyb, trunc=trunc)  # naive reference implementation with N² gate operations (diagonal)
```

> For diagonal coupling, `hybriddynamics` can use either `PartialIF` (site-by-site multiplication, each factor having bond dimension 2) or `XTRGIF`. For off-diagonal coupling, only `XTRGIF` can be used.

## Baths and correlation functions (`ImpurityModelBase`)

```julia
spec = spectrum(f, lb=0, ub=wc)          # define the spectral density J(ω) from a function
spec = Leggett(d=1, ωc=5, α=0.1)         # predefined Leggett spectrum
bath = bosonicbath(spec, β=β)            # bosonic bath (finite temperature β)
corr = correlationfunction(bath, lattice) # correlation function Δ(τ,τ') discretized on the lattice
```

`correlationfunction` automatically calls `Δt` (real time), `Δτ` (imaginary time), or `Δm` (mixed contour) depending on the lattice contour.

Key points:
- zero temperature corresponds to `β=Inf`;
- the integration interval `[lb, ub]` of `spectrum` should cover the main weight of the spectrum (e.g. take `ub=3~5ωc` for an Ohmic spectrum); if $J(\omega)$ diverges at $\omega=0$ or is convolved with a divergent kernel, start from a nonzero lower bound (see the [Practice guide](@ref) for details).

## ADT / PT tensor operations

- `mult(a, b, algmult)` / `mult!(a, b, trunc)`: tensor network multiplication (elementwise product / MPO multiplication);
- `integrate(lattice, args...)` / `integrate(mpsA, mpsB)`: computes the partition function (path-integral sum);
- `apply!(term, mps)`: applies a local operator;
- `canonicalize!`, `leftorth!`, `rightorth!`: orthogonalization (`Orthogonalize` can be specified);
- `bond_dimension(mps)`, `bond_dimensions(mps)`: bond dimension queries;
- `distance(mps1, mps2)` / `distance2`: distance between two tensors (for relative-error validation);
- `randomadt` / `randompt`: random tensors (for testing); `vacuumstate(lattice)`: vacuum state.

## [Observables](@id manual_observables)

There are two measurement paths, both applicable to ADT and PT:

**Path A: operator insertion (arbitrary operators, including off-diagonal and two-point correlations)**

Insert the operator (which can be any matrix) into the system dynamics as a `ContourOperator`, then evaluate everything at once:

```julia
# Single point: ⟨op(t_i)⟩
ct   = ContourOperator(ContourIndex(i), op)
mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)    # real time requires ρ₀; imaginary time does not
mps2 = mult!(mpsK, mpsI, trunc=trunc)
v    = integrate(mps2) / integrate(mps)

# Two-point correlation: ⟨op2(t_i) op1(t_j)⟩
ct   = ContourOperator([ContourIndex(i), ContourIndex(j)], [op2, op1])
mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
# ... likewise: boundarycondition! → mult! → integrate/integrate
```

**Path B: environment caches (batched single-point observables)**

Build the cache once, then measure point by point (already normalized):

```julia
# PT framework (off-diagonal coupling, or off-diagonal observables needed)
cache = environments(lattice, mps, ρ₀=ρimp)
v = expectationvalue(ContourOperator(ContourIndex(i, :+), op), cache)

# ADT framework (diagonal coupling, diagonal observables)
cache = environments(lattice, mpsK, mpsI)
v = expectationvalue(ADTTerm(index(lattice, i, branch=:+), zdiag), cache)
```

Imaginary time / mixed contours are supported as well (PT uses `environments(lattice, mps)`, ADT uses `environments(lattice, mpsK, mpsI)`). The multi-point ADT form `ADTTerm((pos2, pos1), (v2, v1))` can measure diagonal two-point correlations.

Auxiliary functions: `Zvalue(cache)` (partition function), `Zvalue2(cache)`, `TransferMatrix` (transfer matrix),
`correlation(lattice, model, op, mpsI[, ρ0])` (two-point correlation function), `heatcurrents`.

## MPO Hamiltonians (long-range interaction tools)

`src/mpohamiltonian/` provides a standalone tool for constructing MPOs of Hamiltonians with long-range (exponentially decaying) interactions, intended for other one-dimensional quantum many-body problems:

```julia
# SchurMPOTensor: encodes [local terms + exponentially decaying long-range terms] into a compact MPO site tensor
h = SchurMPOTensor(h1, h2s)    # h2s is a list of ExponentialDecayTerm / GenericDecayTerm / PowerlawDecayTerm
mpo = MPOHamiltonian([h for _ in 1:L])
tensors = tompotensors(mpo)              # convert to dense MPO site tensors
tensors2 = timeevompo(tensors, dt, WII())   # time evolution (WI / WII / ComplexStepper / FirstOrderStepper)
```

## Hyperparameters and sources of error

The method has four classes of controllable error sources (Sec. IV of the paper):

| Hyperparameter | Meaning | Suggested default | Error source |
|---|---|---|---|
| `δt` / `δτ` | Contour discretization step | determined by convergence checks | QuAPI first-order Trotter decomposition error, overall $O(t\delta t)$ |
| `χ` (`trunc.D`) | MPO/MPS bond dimension truncation | 30–200 in the paper's examples | MPO compression (truncation) error |
| `n` (`algexpan`) | Number of terms in the exponential expansion of the hybridization function | 20 by default in the paper | Prony approximation error |
| `m` (`alg.k`) | Number of XTRG steps (time step $1/2^m$) | 7 by default in the paper | XTRG thermal-state construction error (exponentially convergent) |
| `d` (lattice) | Local Hilbert space truncation of the bosonic impurity | depends on temperature | local truncation error |

Key points (all verified in the paper with examples such as the JC model, two single-mode baths, non-interacting bosonic modes, and a sub-Ohmic bath):

- errors saturate quickly as `χ`, `m`, and `n` increase; the error decreases monotonically as `δt` shrinks;
- strong coupling is generally harder to converge than weak coupling (requires larger `χ`);
- JC-type (number-conserving) coupling generates far less entanglement than Rabi-type coupling and usually converges faster;
- the zero-temperature case ($\beta=\infty$) is generally easier to converge than the finite-temperature case;
- the computational cost is dominated by MPO–MPO multiplications in XTRG, scaling theoretically as $O(N\chi^4 d^3)$.

## Code structure

```
src/
├── TEMPO.jl                    # module definition and exported symbols
├── tensorops/                  # basic tools for truncation, tensor operations, orthogonalization, etc. (algorithms.jl contains DMRG multiplication / MPS algorithms)
├── defaults.jl                 # default hyperparameters
├── mpohamiltonian/             # MPO Hamiltonians (SchurMPO, long-range terms, time-evolution steppers)
├── contourindices.jl           # ContourIndex, branch
├── adt/  pt/                   # augmented density tensor (MPS) and process tensor (MPO) data structures and operations
├── conversions.jl              # PT ↔ ADT conversion helpers (copy tensors etc.)
├── adtterms.jl / fockterms.jl  # local terms such as ADTTerm / FockTerm / ProdFockTerm
├── adtlattices/ ptlattices/    # real / imaginary / mixed contour lattice definitions (Fock ordering)
├── contouroperators.jl         # ContourOperator (operators on the PT)
├── correlationfunction.jl      # discretization of bath correlation functions onto the lattice
├── influencefunctional/        # Feynman-Vernon IF: PartialIF / translationally invariant IF (both ADT and PT)
├── tdinfluencefunctional/      # IF for time-dependent couplings and time-dependent hybridization styles
├── boundarycondition.jl        # initial state / boundary conditions
├── models/                     # unitary (ImpurityHamiltonian) and dissipative (ImpurityLindbladian) dynamics
└── observables/                # environment caches, expectation values, transfer matrices, correlation functions, heat currents

docs/tutorials/                  # tutorials: spinboson and three paper-reproduction notebooks
benchmark/                       # benchmark cases: spin–boson, independent bosonic modes, real time, etc.
test/                            # test suite (including cross-checks against exact diagonalization)
```

## Correspondence with the literature

| Paper | Package implementation |
|---|---|
| Diagonal + commuting coupling (standard TEMPO) | `ADTLattice` + `AdditiveHyb` + `PartialIF`/`XTRGIF` |
| Diagonal + non-commuting coupling (multiple baths) | IFs of multiple `AdditiveHyb` multiplied together (`mult!`) |
| Off-diagonal coupling (conjugate pair, the core generalization) | `PTLattice` + `NonDiagonalHyb` + `XTRGIF` |
| Real-time Keldysh contour | `contour=:real` |
| Imaginary-time contour | `contour=:imag` |
| L-shaped Kadanoff–Baym contour | `contour=:mixed` (`MixedPTLattice`/`MixedADTLattice`) |
| QuAPI discretization (Appendix C) | `Δt`/`Δτ`/`Δm` in `correlationfunction(bath, lattice)` |
| Exponential expansion of the hybridization function (Appendix D) | `OverDeterminedProny`/`DeterminedProny` (`algexpan`) |
| XTRG construction of the effective thermal state (Fig. 3 of the paper) | `XTRGIF(k=..., fast=...)` |
| Absorption of the system Hamiltonian into the PT (Fig. 2b of the paper) | `sysdynamics(lattice, model)` + `mult!` |
| Observable computation (Fig. 1d, 1e of the paper) | `environments` + `expectationvalue` |
| JC spin–boson model | `docs/tutorials/spinboson/jctype.jl` |
| Standard spin–boson model | `docs/tutorials/spinboson/rabitype.jl`, `benchmark/sb.jl` |
| Non-interacting / interacting bosonic impurity | `benchmark/independentbosons.jl`, `bosonicimpurity.jl` |
