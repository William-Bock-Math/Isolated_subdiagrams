# ==============================================================================
#  Isolated subdiagrams of the maximal affine subdiagrams, via k-odd graphs.
#
#  Needs the packages Oscar and Combinatorics.
# ==============================================================================

using Oscar

include("k_odd_fast.jl")
include("Diagrams.jl")

"""
    injective_on_components(D, P_aff, k, c_type; P = nothing, comps_D = …)

Compares the k-odd graph `P^k` of the maximal affine subdiagram `P_aff` with the
k-odd graph `D^k` of the whole diagram, both for inc-subdiagrams of type
`c_type`.  `P^k` is an induced subgraph of `D^k`, so each component of `P^k`
lies in a unique component of `D^k`; it is *isolated* if no other component of
`P^k` lies there.

Returns one named tuple `(component, host, count)` per isolated component: the
component of `P^k` (if `P` is given, only its subdiagrams contained in `P`),
the component of `D^k` containing it, and `length(component)`.

`comps_D` defaults to `odd_components(D, k; c_type)`; pass it in to reuse it
for several `P_aff`.  `D` is `diagram_data(labels)` or the `labels` dict.
"""
function injective_on_components(D::NamedTuple, P_aff, k::Int, c_type; P = nothing,
                                 comps_D = odd_components(D, k; c_type))
    comps_P = odd_components(D, k; c_type, subdiag = P_aff)

    host  = Dict(I => h for (h, comp) in enumerate(comps_D) for I in comp)
    hosts = map(comps_P) do comp
        h = host[first(comp)]
        @assert all(host[I] == h for I in comp)    # P^k is a subgraph of D^k
        h
    end

    n_inside = zeros(Int, length(comps_D))         # components of P^k per component of D^k
    for h in hosts
        n_inside[h] += 1
    end

    isolated = [(comp, h) for (comp, h) in zip(comps_P, hosts) if n_inside[h] == 1]
    return map(isolated) do (comp, h)
        sub = P === nothing ? comp : [I for I in comp if I ⊆ P]
        (component = sub, host = comps_D[h], count = length(sub))
    end
end
injective_on_components(labels::AbstractDict, args...; kwargs...) =
    injective_on_components(diagram_data(labels), args...; kwargs...)

"""
    IsolatedCounts

The result of `isolated_subdiagrams`, displayed as a table.  `r.counts[k, j]`
(or `r[k, j]`) is the number for rank `k` and the maximal affine subdiagram
`affines[j]`, whose column is labelled `r.names[j]`.
"""
struct IsolatedCounts
    counts::Matrix{Int}
    names::Vector{String}
end

Base.getindex(r::IsolatedCounts, k::Integer, j::Integer) = r.counts[k, j]
Base.show(io::IO, r::IsolatedCounts) = print(io, "IsolatedCounts(", join(size(r.counts), " × "), ")")

function Base.show(io::IO, ::MIME"text/plain", r::IsolatedCounts)
    cell(x) = iszero(x) ? "·" : string(x)
    wk = max(1, ndigits(size(r.counts, 1)))
    w  = [max(textwidth(name), maximum(textwidth ∘ cell, col; init = 1))
          for (name, col) in zip(r.names, eachcol(r.counts))]
    println(io, "Isolated subdiagrams by rank k and maximal affine subdiagram (· = 0):")
    println(io)
    println(io, " ", lpad("k", wk), " │ ", join(lpad.(r.names, w), "   "))
    print(io, "─"^(wk + 2), "┼", "─"^(sum(w; init = 0) + 3 * length(w) - 1))
    for (k, row) in enumerate(eachrow(r.counts))
        print(io, "\n ", lpad(k, wk), " │ ", join(lpad.(cell.(row), w), "   "))
    end
end

"""
    isolated_subdiagrams(labels, dimension, affines, affines_par; names = nothing)

For every rank `k = 1, …, dimension` and every maximal affine subdiagram
`affines[j]`, the total `count` of its isolated components (see
`injective_on_components`, with `P = affines_par[j]`), summed over all inc types
of rank `k` that occur in the diagram.

Returns an `IsolatedCounts`, which displays as a table with one row per rank and
one column per affine subdiagram; the columns are labelled by the types of
`affines_par`, or by `names` if given.  In a script, show it with `display`.
"""
function isolated_subdiagrams(labels, dimension::Int, affines, affines_par; names = nothing)
    length(affines) == length(affines_par) ||
        throw(DimensionMismatch("affines and affines_par must have the same length"))
    D = diagram_data(labels)
    if names === nothing
        names = map(enumerate(affines_par)) do (j, P)
            t = coxeter_type(D.M, P)
            t === nothing ? "#$j" : type_string(t)
        end
    end
    length(names) == length(affines) ||
        throw(DimensionMismatch("need one name per affine subdiagram"))

    counts = zeros(Int, dimension, length(affines))
    for k in 1:dimension
        for t in inc_types_D(D, k)
            comps_D = odd_components(D, k; c_type = t)    # the same for every affine
            for (j, (P_aff, P)) in enumerate(zip(affines, affines_par))
                rows = injective_on_components(D, P_aff, k, t; P, comps_D)
                counts[k, j] += sum((r.count for r in rows); init = 0)
            end
        end
        delete!(D.spherical, k)                           # rank k is done, free the cache
    end
    return IsolatedCounts(counts, string.(names))
end
