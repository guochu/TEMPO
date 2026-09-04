# TEMPO.jl

This package is a Julia implementation of the **Time-Evolving Matrix Product Operator (TEMPO)** method, with the algorithmic theory based on the publication:

> C. Guo, W. Wu, X. Xu, T. Jiang, P.-X. Chen, and R. Chen,
> *Time-evolving matrix product operators for off-diagonal system-bath coupling*,
> **Phys. Rev. B 114, 125413 (2026)**.

Unlike the original TEMPO, which only supports diagonal system-bath coupling, this implementation is built on the **process tensor (PT)** framework and generalizes TEMPO to the more general case of **off-diagonal system-bath coupling** (the system couples to the bath through a pair of conjugate non-Hermitian operators $\hat{A}^\dagger, \hat{A}$), while unifying:

- standard TEMPO (diagonal, commuting coupling, `ADT` + partial influence functional);
- diagonal but non-commuting multi-bath coupling;
- off-diagonal (conjugate-pair) coupling (`PT` + translationally invariant influence functionals);
- evolution on real-time, imaginary-time, and mixed (Kadanoff–Baym) contours;
- time-dependent system-bath coupling.

## Documentation

| Page | Content |
|---|---|
| [Quickstart](@ref) | Complete workflows for three classes of typical problems (can be copied and run as-is) |
| [Manual](@ref) | Core components, hyperparameters and error sources, code structure, correspondence with the literature |
| [Practice guide](@ref) | Practical experience: path selection, convergence checks, common pitfalls and diagnostics (**recommended reading before any general computation task**) |
| [Internals](@ref) | Internal data structures, algorithm flow, and tensor conventions |
| [API reference](@ref) | Item-by-item documentation of all exported symbols |

Notebooks reproducing the accompanying paper can be found in `docs/tutorials/` (strathearn2018, otterpohl2025, guo2026, spinboson).

## Installation

This package depends on the following packages (see `Project.toml`):

| Package | Purpose |
|---|---|
| `ImpurityModelBase` | Defines basic objects such as spectral densities, baths (`bosonicbath`), and bosonic operators |
| `QuAPI` | Provides basic contour types such as `ContourIndex`, `branch`, and `index` |
| `ExpExp` | Exponential expansion of hybridization functions (Prony method) |
| `TensorOperations` | Tensor network contractions |
| `Polynomials`, `LinearAlgebra`, `Statistics`, `TupleTools`, `Logging` | General-purpose numerical and linear-algebra utilities |

After activating the project in Julia, it can be used as follows:

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra
```

!!! note
    Functions such as `spectrum`, `bosonicbath`, `bosonaoperator`, `bosondensityoperator`, and `Leggett` are defined in `ImpurityModelBase` and must be loaded together with `using`.

## Core concepts

### The quantum impurity problem (QIP)

Consider an "impurity" system $\hat{H}_S$ linearly coupled to a non-interacting bosonic bath:

```math
\hat{H} = \hat{H}_S + \hat{H}_{\text{int}},
\qquad
\hat{H}_{\text{int}} = \hat{H}_{\text{hyb}} + \hat{H}_B .
```

- **Diagonal coupling** (original TEMPO/QuAPI):

  ```math
  \hat{H}_{\text{hyb}} = \sum_{l,k} \hat{A}_l\, (V_{l,k} \hat{b}^\dagger_{l,k} + \mathrm{H.c.}),
  ```
  where $\hat{A}_l$ is a Hermitian operator, and in the coupling term $\hat{A}_l$ only appears paired with linear combinations of $\hat{b}^\dagger+\hat{b}$.

- **Off-diagonal coupling** (the generalization introduced in this work):

  ```math
  \hat{H}_{\text{hyb}} = \sum_{l,k} \left( V_{l,k} \hat{A}_l \hat{b}^\dagger_{l,k} + \mathrm{H.c.} \right),
  ```
  where $\hat{A}_l$ can be a non-Hermitian operator (e.g., the Jaynes–Cummings type $\hat{A}=\hat{\sigma}_-$).

Off-diagonal coupling cannot be recombined into a "diagonal but non-commuting" case and therefore requires a new framework.

### The Feynman–Vernon influence functional (IF)

The key starting point of TEMPO-type methods is the Feynman–Vernon influence functional obtained after tracing out the bath. For off-diagonal coupling, it takes the operator path form on the Keldysh contour $C$:

```math
\mathcal{I}[\hat{A}^\dagger, \hat{A}]
= \mathcal{T}_C \exp\left[ -\int_C \mathrm{d}\tau \int_C \mathrm{d}\tau'\,
    \hat{A}^\dagger(\tau)\, \Delta(\tau,\tau')\, \hat{A}(\tau') \right],
```

where the **hybridization function** is determined by the spectral density:

```math
\Delta(\tau,\tau') = i \int \mathrm{d}\omega\, J(\omega)\, D_\omega(\tau,\tau'),
\qquad
J(\omega) = \sum_k |V_k|^2 \delta(\omega - \omega_k).
```

For diagonal coupling, this IF corresponds to a classical partition function (and can be represented as an MPS / ADT); for off-diagonal coupling, it corresponds to the thermal state $\mathrm{e}^{-\hat{H}_{\text{eff}}}$ of an effective quantum many-body Hamiltonian and must be represented as an **MPO (i.e., a PT)**.

### Two tensor networks: ADT and PT

| Object | Full name | Representation | Applicable coupling |
|---|---|---|---|
| `ADT` | Augmented density tensor | MPS | Diagonal coupling (original TEMPO) |
| `ProcessTensor` (`PT`) | Process tensor | MPO | Diagonal non-commuting and off-diagonal coupling (extension of this work) |

A PT can be systematically converted into an ADT by inserting 3D copy tensors between neighboring sites (Fig. 4 of the paper).

### Time contours

This package supports three time contours, selected via the `contour` keyword:

- `contour=:real` (equivalent to `:Keldysh`): real-time evolution, with the initial system state $\hat{\rho}_S \otimes \hat{\rho}_B$;
- `contour=:imag`: imaginary-time evolution ($0\to\beta$), corresponding to the finite-temperature partition function and Matsubara correlation functions;
- `contour=:mixed` (equivalent to `:Kadanoff`): the L-shaped Kadanoff–Baym contour (a mix of imaginary and real time).
