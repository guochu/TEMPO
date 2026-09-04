# TEMPO.jl

**TEMPO.jl** is a Julia implementation of the [time-evolving matrix product operator (TEMPO)](https://link.aps.org/doi/10.1103/PhysRevLett.118.077601) method for the accurate simulation of **open quantum impurity problems** (a quantum system linearly coupled to a continuous bosonic bath characterized by a spectral density). The algorithmic theory is based on:

> C. Guo, W. Wu, X. Xu, T. Jiang, P.-X. Chen, and R. Chen,
> *Time-evolving matrix product operators for off-diagonal system-bath coupling*,
> **Phys. Rev. B 114, 125413 (2026)**.

In contrast to the original TEMPO, which only supports diagonal system-bath coupling, this implementation is built on the **process tensor (PT)** framework and generalizes TEMPO to general **off-diagonal system-bath coupling** (the system couples to the bath through a conjugate pair of operators $\hat{A}^\dagger,\hat{A}$), within a unified framework supporting:

- standard TEMPO (diagonal, commuting coupling, `ADT`/MPS representation + partial influence functional);
- diagonal but non-commuting multi-bath coupling;
- off-diagonal (conjugate-pair) coupling (`ProcessTensor`/MPO representation + XTRG-style translationally invariant influence functional);
- real-time (Keldysh), imaginary-time (Matsubara), and mixed (Kadanoff–Baym) contours;
- dissipative impurities (Lindblad-type system dynamics, `ImpurityLindbladian`).

## Installation

Julia ≥ 1.10 is required. The package depends on three companion packages, all open-sourced on [GitHub (guochu)](https://github.com/guochu): `ImpurityModelBase` (basic objects such as spectral densities, baths, and bosonic operators), `QuAPI` (contour basic types such as `ContourIndex`), and `ExpExp` (Prony exponential expansion of hybridization functions).

```bash
git clone <this repo>            # place it in the same parent directory as ImpurityModelBase and QuAPI
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Usage:

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra
```

See [Project.toml](Project.toml) for the full list of dependencies.

## Quick start

A spin-1/2 impurity coupled diagonally via $\hat{\sigma}_z$ to a sub-Ohmic bath; computing the real-time evolution of $\langle\sigma_z(t)\rangle$ (the code can be copied and run as is):

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra

# 1. truncation scheme and real-time lattice
trunc   = truncdimcutoff(D=32, ϵ=1.0e-12)
lattice = ADTLattice(N=100, δt=0.05, contour=:real)

# 2. system: diagonal σ_z coupling; bath: sub-Ohmic spectrum
z    = Matrix{ComplexF64}([-1 0; 0 1])
hyb  = AdditiveHyb([z[1,1], z[2,2]])                    # pass the diagonal entries of the coupling operator
spec = spectrum(w -> 2π*0.1*w^0.5*5^(1-0.5), lb=0, ub=5) # J(ω) = 2π α ω^s ωc^(1-s)
bath = bosonicbath(spec, β=1.0)
corr = correlationfunction(bath, lattice)

# 3. build the influence functional (IF, an MPS)
mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)

# 4. bare system dynamics + initial-state boundary condition
model = ImpurityHamiltonian(0.2 .* Matrix{ComplexF64}([0 1; 1 0]))
ρ₀    = Matrix{ComplexF64}([1 0; 0 0])                  # σ_z-polarized initial state
mpsK  = boundarycondition!(sysdynamics(lattice, model, trunc=trunc), lattice, ρ₀=ρ₀)

# 5. environment cache and observables ⟨σ_z(t_i)⟩
cache = environments(lattice, mpsK, mpsI)
sz = [expectationvalue(ADTTerm(index(lattice, i, branch=:+), [z[1,1], z[2,2]]), cache) for i in 1:100]
```

Every computation follows the same template:

```text
1. trunc   = truncation scheme
2. lattice = ADTLattice / PTLattice (contour + number of steps + step size)
3. hyb     = AdditiveHyb / NonDiagonalHyb / ... (coupling operators)
4. spec/bath = spectrum + bosonicbath (bath spectrum and temperature)
5. corr    = correlationfunction(bath, lattice)
6. mpsI    = hybriddynamics(lattice, corr, hyb[, alg])   # influence functional
7. mpsK    = sysdynamics(lattice, model) + boundarycondition!  # system dynamics + initial state
8. cache   = environments(...); expectationvalue(...)    # observables
```

For the other typical workflows — off-diagonal (Jaynes–Cummings-type) coupling, imaginary-time evolution, and time-dependent coupling — see the [quickstart](docs/src/quickstart.md).

## Documentation

The full documentation is built with Documenter.jl: `julia --project=docs docs/make.jl`.

| Document | Content |
|---|---|
| [Introduction and core concepts](docs/src/index.md) | Method background: quantum impurity problem, influence functional, ADT vs. PT, time contours |
| [Quickstart](docs/src/quickstart.md) | Complete runnable workflows for four typical problems (diagonal/off-diagonal coupling, imaginary time, time-dependent coupling) |
| [Manual](docs/src/manual.md) | Core components, hyperparameters and error sources, code structure, correspondence with the literature |
| [Practice guide](docs/src/practice.md) | Path selection, convergence checks, common pitfalls and diagnostics (recommended reading before general-purpose computations) |
| [Internals](docs/src/internals.md) | Internal data structures, algorithmic workflow, and tensor conventions |
| [API reference](docs/src/api.md) | Itemized documentation of all exported symbols |
| [Interface changes](changes.md) | Recent API changes (tensorops restructuring, `tsvd`/`leftorth`/`rightorth`, etc.) |

## Code structure

```
src/
├── TEMPO.jl                    # module definition and exported symbols
├── tensorops/                  # basic tools: truncation schemes, tensor factorizations (tsvd/leftorth/rightorth), tensor operations
├── algorithms.jl               # MPS/DMRG algorithms (SVDCompression etc.)
├── defaults.jl                 # default hyperparameters
├── mpohamiltonian/             # MPO Hamiltonians (SchurMPO, long-range decay terms, time-evolution steppers)
├── adt/  pt/                   # ADT (MPS) and process tensor (MPO) data structures and operations
├── adtlattices/  ptlattices/   # real/imaginary/mixed-time lattices and index mappings
├── correlationfunction.jl      # bath correlation functions
├── influencefunctional/        # influence functionals: PartialIF / XTRGIF / TDVPIF / PT paths
├── tdinfluencefunctional/      # time-dependent system-bath coupling
├── boundarycondition.jl        # initial states / boundary conditions
├── models/                     # system Hamiltonians / Lindblad operators
└── observables/                # environment caches, expectation values, partition functions, transfer matrices
```

## Testing

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Covers all core paths (lattices, influence functional construction, time-dependent coupling, model benchmarks, observables, etc.).

## Tutorials and paper reproductions

`docs/tutorials/` contains reproduction notebooks for several papers (with model descriptions, parameter documentation, and MPS-IF caches):

- `strathearn2018` — original TEMPO (PRL 120, 060601 (2018))
- `otterpohl2025` — PT representation for off-diagonal coupling
- `guo2026` — this paper (Phys. Rev. B 114, 125413 (2026))
- `spinboson` — Rabi-type / JC-type spin-boson models

Generation and execution:

```bash
julia --project=docs docs/tutorials/gen_notebooks.jl   # generate .ipynb files
julia --project=docs docs/tutorials/run_notebook.jl    # run without a Jupyter environment
```

## Related project: GTEMPO

For **fermionic** impurity problems, check out our sister package [GTEMPO](https://github.com/guochu/GTEMPO), which implements the Grassmann time-evolving matrix product operator (GTEMPO) method on the Keldysh, imaginary-time, and Kadanoff contours, together with fermionic QuAPI, translationally invariant IF algorithms, current measurements, and more.

## Todo list

- [ ] Support for time-translational invariance (infinite MPS, see e.g. PRL 132, 200403 (2024))
- [ ] Partial integration technique for large impurities (Phys. Rev. B 112, 155115 (2025))
- [ ] Support for the continuous MPS technique (PRB 110, 045104 (2024))
- [ ] Time-dependent system-reservoir coupling
- [ ] Higher-order QuAPI

## Citation

If this package is helpful to your research, please cite:

```bibtex
@article{GuoChen2026,
  title   = {Time-evolving matrix product operators for off-diagonal system-bath coupling},
  author  = {Guo, C. and Wu, W. and Xu, X. and Jiang, T. and Chen, P.-X. and Chen, R.},
  journal = {Physical Review B},
  volume  = {114},
  pages   = {125413},
  year    = {2026},
  doi     = {10.1103/PhysRevB.114.125413}
}
```
