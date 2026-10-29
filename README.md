# Isolated subdiagrams and k-odd graphs of Coxeter diagrams

Julia code, built on [Oscar](https://www.oscar-system.org), for computing the
k-odd graphs of a Coxeter diagram (in the sense of Michael, Santos Rego, Schwer
and Varghese) and counting isolated subdiagrams in its maximal affine
subdiagrams. Edge labels can be arbitrary, including ∞.

<!-- TODO: link the paper, e.g. [arXiv:XXXX.XXXXX](https://arxiv.org/abs/XXXX.XXXXX) -->

## Installation

You need [Julia](https://julialang.org/downloads/) 1.10 or newer. Then run

```sh
git clone https://github.com/USERNAME/isolated-subdiagrams.git
cd isolated-subdiagrams
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This installs Oscar and Combinatorics into an environment that belongs to this
folder. The first time, installing and precompiling Oscar takes several minutes.
Without git, use *Code → Download ZIP* on this page instead of `git clone`.

## Quick start

Start Julia in the repository folder with `julia --project=.`, then:

```julia
julia> include("isolated_subdiagrams.jl")

julia> isolated_subdiagrams(labels_HE8, 9, affines_HE8, affines_par_HE8)
Isolated subdiagrams by rank k and maximal affine subdiagram (· = 0):

 k │ E₈
───┼────
 1 │  8
 2 │ 21
 3 │ 21
 4 │  1
 5 │  2
 6 │  1
 7 │  1
 8 │  1
 9 │  ·
```

[`examples.jl`](examples.jl) has more sample computations. Run them all with
`julia --project=. examples.jl`.

## Your own diagrams

A diagram on the vertices 1, …, n is a dictionary of Coxeter labels
`(i, j) => m`, where `m` is the order of sᵢsⱼ: an integer ≥ 3, or `Inf`
for infinite order. Pairs that are not listed commute.

```julia
labels = Dict((1, 2) => 3, (2, 3) => 4, (3, 4) => 4)    # the hyperbolic group [3, 4, 4]
isolated_subdiagrams(labels, 3, [[2, 3, 4]], [[2, 3]])
```

The arguments of `isolated_subdiagrams(labels, dimension, affines, affines_par)`
are the diagram, the largest rank `k` to compute, the maximal affine
subdiagrams (as vertex lists), and for each of them a maximal spherical
subdiagram inside it. The columns of the table are labelled by the types of
`affines_par`; pass `names = [...]` to choose your own.
[`Diagrams.jl`](Diagrams.jl) contains the diagrams used in the examples.

## What is computed

* **Spherical subdiagrams** are recognised, and their types determined, by the
  classification of finite Coxeter groups: Aₙ, Bₙ, Dₙ, E₆, E₇, E₈, F₄, G₂, H₃,
  H₄, I₂(m).
* **Inc types** are the irreducible types whose longest element is central:
  A₁, Bₙ, Dₙ (n even), E₇, E₈, F₄, G₂, H₃, H₄ and I₂(m) (m even). An *inc
  subdiagram* `I` is a spherical subdiagram all of whose components have inc
  type, and `X_I` is the union of its components that are not of type A₁.
* **The k-odd graph** `D^k` has the inc subdiagrams of rank k as vertices;
  `I` and `J` are adjacent if `X_I = X_J` and `J = (I ∖ {s}) ∪ {t}` for an edge
  `{s, t}` with an odd label `m_st`.
* **`isolated_subdiagrams`**: take a maximal affine subdiagram `P_aff` with
  spherical part `P`, and an inc type. The k-odd graph of `P_aff` is a subgraph
  of `D^k`, and one of its components is *isolated* if it is the only one inside
  its component of `D^k`. The table entry for rank k is the number of
  subdiagrams contained in `P` that lie in isolated components, summed over all
  inc types of rank k.

## Files

| File | Contents |
|---|---|
| [`isolated_subdiagrams.jl`](isolated_subdiagrams.jl) | entry point (`include` this): isolated components and the result table |
| [`k_odd_fast.jl`](k_odd_fast.jl) | Coxeter types, inc subdiagrams, k-odd graphs |
| [`Diagrams.jl`](Diagrams.jl) | the example diagrams (K3, HE8, …) |
| [`examples.jl`](examples.jl) | sample computations |

The main functions are `isolated_subdiagrams`, `injective_on_components`,
`k_odd`, `odd_components`, `inc_subdiagrams`, `inc_types_D`, `coxeter_type`,
`type_string`, `diagram_data` and `coxeter_matrix_from_labels`. All of them
have docstrings: type for example `?k_odd` in the Julia REPL.

## Performance

Spherical subdiagrams are found by checking every vertex subset of each rank,
so the running time roughly doubles with each additional vertex. The 24-vertex
K3 diagram is the largest example here.

## License

See [LICENSE](LICENSE).
