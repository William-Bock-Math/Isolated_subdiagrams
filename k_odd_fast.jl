# ==============================================================================
#  The k-odd graph D^k  (Michael – Santos Rego – Schwer – Varghese)
#
#  A Coxeter diagram on the vertices 1, …, n is given by a dictionary `labels`
#  of Coxeter labels (i, j) => m_ij, the order of s_i s_j: an integer ≥ 3, or
#  Inf for infinite order.  Pairs that are not listed commute (m_ij = 2;
#  listing them with 2 is also fine).  The order of i and j in a key does not
#  matter, and n is the largest vertex index that occurs.  Example:
#
#      labels = Dict((1, 2) => 3, (2, 3) => 5, (3, 4) => Inf)
#
#  Originally AI-generated, cleaned up.
# ==============================================================================

using Oscar
import Combinatorics: combinations

# ------------------------------------------------------------------------------
# 1. Coxeter matrices and types
# ------------------------------------------------------------------------------

"""
    coxeter_matrix_from_labels(labels; n = largest vertex index)

The Coxeter matrix `M` as a `Matrix{Int}`: `M[i, i] = 1`, `M[i, j] = m_ij`, and
`M[i, j] = 0` for infinite order (the convention of Oscar's `coxeter_matrix`).
"""
function coxeter_matrix_from_labels(labels::AbstractDict; n::Int = maximum(maximum, keys(labels)))
    M = fill(2, n, n)
    for i in 1:n
        M[i, i] = 1
    end
    for ((i, j), m) in labels
        i != j || throw(ArgumentError("label on the loop ($i, $i)"))
        if isinf(m) && m > 0
            mij = 0
        elseif isinteger(m) && m >= 2
            mij = Int(m)
        else
            hint = m == 1 ? " (old convention? see convert_old_labels)" : ""
            throw(ArgumentError("label $m on ($i, $j): must be an integer ≥ 2 or Inf$hint"))
        end
        M[i, j] == 2 || M[i, j] == mij ||
            throw(ArgumentError("conflicting labels for ($i, $j) and ($j, $i)"))
        M[i, j] = M[j, i] = mij
    end
    return M
end

"""
    convert_old_labels(old)

Translates labels in the old convention (`1` = order 3, `2` = infinite order) to
Coxeter labels.
"""
convert_old_labels(old::AbstractDict) =
    Dict(e => (w == 1 ? 3 : Inf) for (e, w) in old if !iszero(w))

"""
    coxeter_type(M, I)

The type of the subdiagram on the vertices `I` of the diagram with Coxeter
matrix `M`: the sorted vector of its irreducible types, or `nothing` if it
generates an infinite group.  The irreducible types are `(:A, n)`, `(:B, n)`,
`(:D, n)`, `(:E, n)`, `(:F, 4)`, `(:G, 2)`, `(:H, n)`, and `(:I, m)` for the
dihedral group `I_2(m)` with `m = 5` or `m ≥ 7`.
"""
function coxeter_type(M::AbstractMatrix{<:Integer}, I::AbstractVector{<:Integer})
    k = length(I)
    k == 0 && return Tuple{Symbol,Int}[]         # the empty diagram
    edges = 0
    for a in 1:k, b in a+1:k                     # cheap rejections first
        m = M[I[a], I[b]]
        m == 2 && continue
        m == 0 && return nothing                 # infinite order
        edges += 1
    end
    edges < k || return nothing                  # finite diagrams are forests

    nbrs  = [[b for b in 1:k if b != a && M[I[a], I[b]] != 2] for a in 1:k]
    label = (a, b) -> M[I[a], I[b]]
    types = Tuple{Symbol,Int}[]
    seen  = falses(k)
    for a in 1:k
        seen[a] && continue
        comp = [a]                               # the component of a; grows while
        seen[a] = true                           # the loop below runs over it
        for v in comp, b in nbrs[v]
            seen[b] || (seen[b] = true; push!(comp, b))
        end
        t = irreducible_type(comp, nbrs, label)
        t === nothing && return nothing
        push!(types, t)
    end
    return sort!(types)
end

# Classification of the connected Coxeter diagrams of finite type.  `comp` is a
# connected component, `nbrs` its adjacency lists, `label(a, b)` the edge labels.
function irreducible_type(comp, nbrs, label)
    c = length(comp)
    c == 1 && return (:A, 1)
    deg(a) = length(nbrs[a])
    sum(deg, comp) == 2(c - 1) || return nothing     # contains a cycle
    if c == 2
        m = label(comp[1], comp[2])
        return m == 3 ? (:A, 2) : m == 4 ? (:B, 2) : m == 6 ? (:G, 2) : (:I, m)
    end

    heavy  = [(a, b) for a in comp for b in nbrs[a] if a < b && label(a, b) > 3]
    branch = [a for a in comp if deg(a) > 2]
    if isempty(heavy)                                # simply laced: A_n, D_n, E_n
        isempty(branch) && return (:A, c)
        length(branch) == 1 && deg(branch[1]) == 3 || return nothing
        arms = sort!([arm_length(nbrs, branch[1], a) for a in nbrs[branch[1]]])
        arms[1] == 1 || return nothing
        arms[2] == 1 && return (:D, c)
        arms[2] == 2 && arms[3] <= 4 && return (:E, c)
        return nothing
    end

    # a path with one heavy edge: B_n, F_4, H_3, H_4
    length(heavy) == 1 && isempty(branch) || return nothing
    a, b = heavy[1]
    m = label(a, b)
    at_end = deg(a) == 1 || deg(b) == 1
    m == 4 && at_end && return (:B, c)
    m == 4 && c == 4 && return (:F, 4)
    m == 5 && at_end && c <= 4 && return (:H, c)
    return nothing
end

# number of vertices on the arm that starts at v, coming from the branch vertex
function arm_length(nbrs, from, v)
    len = 1
    while length(nbrs[v]) == 2
        from, v = v, (nbrs[v][1] == from ? nbrs[v][2] : nbrs[v][1])
        len += 1
    end
    return len
end

"""
    type_string(t)

A type such as `[(:A, 1), (:E, 8), (:E, 8)]` written as `"A₁E₈²"`.
"""
function type_string(t)
    isempty(t) && return "∅"
    sub(n) = join(('₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉')[d + 1] for d in reverse(digits(n)))
    sup(n) = join(('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')[d + 1] for d in reverse(digits(n)))
    return join(map(unique(t)) do (X, n)
        name = X === :I ? "I₂($n)" : string(X, sub(n))
        e = count(==((X, n)), t)
        e == 1 ? name : name * sup(e)
    end)
end

"""
    elliptic_subdiagrams(M, n; base = nothing)

All vertex sets `I` of size `n` that generate a finite group, i.e. the spherical
subdiagrams of rank `n`, as sorted vectors in lexicographic order; `M` is the
Coxeter matrix.  With `base = v`, only those containing `v`.
"""
function elliptic_subdiagrams(M::AbstractMatrix{<:Integer}, n::Int; base = nothing)
    n > 0 || return [Int[]]
    return [I for I in combinations(1:size(M, 1), n)
            if (base === nothing || base in I) && coxeter_type(M, I) !== nothing]
end

# ------------------------------------------------------------------------------
# 2. The diagram
# ------------------------------------------------------------------------------

"""
    diagram_data(labels)

Everything that depends only on the diagram: its Coxeter matrix `M`, the
neighbour lists `nbrs` (all edges) and `oddnbrs` (edges with odd `m_ij`), and a
cache of the spherical subdiagrams of each rank.

Build it once and pass it to the functions below.  They also accept `labels`
directly, but then nothing is reused between calls.
"""
function diagram_data(labels::AbstractDict)
    M = coxeter_matrix_from_labels(labels)
    n = size(M, 1)
    nbrs    = [[j for j in 1:n if j != i && M[i, j] != 2] for i in 1:n]
    oddnbrs = [[j for j in 1:n if j != i && isodd(M[i, j])] for i in 1:n]
    spherical = Dict{Int,Vector{Pair{Vector{Int},Vector{Tuple{Symbol,Int}}}}}()
    return (M = M, nbrs = nbrs, oddnbrs = oddnbrs, spherical = spherical)
end

"""
    spherical_subdiagrams(D, n; subdiag = nothing)

The spherical subdiagrams `I` of rank `n`, each paired with its type as
`I => type`; with `subdiag`, only those with `I ⊆ subdiag`.  This is a brute
force search over all `binomial(number of vertices, n)` subsets, so it is
computed once per rank and cached in `D`.
"""
function spherical_subdiagrams(D::NamedTuple, n::Int; subdiag = nothing)
    all_n = get!(D.spherical, n) do
        typed = Pair{Vector{Int},Vector{Tuple{Symbol,Int}}}[]
        for I in combinations(1:size(D.M, 1), n)
            t = coxeter_type(D.M, I)
            t === nothing || push!(typed, I => t)
        end
        typed
    end
    subdiag === nothing && return all_n
    S = Set(subdiag)
    return filter(p -> all(in(S), first(p)), all_n)
end

# ------------------------------------------------------------------------------
# 3. inc types and inc subdiagrams
# ------------------------------------------------------------------------------

"""
    is_inc_type((X, n))

Whether the longest element of the irreducible finite type `X_n` is central, i.e.
acts as `-1`: true for `A_1`, `B_n`, `D_n` with `n` even, `E_7`, `E_8`, `F_4`,
`G_2`, `H_3`, `H_4` and `I_2(m)` with `m` even.  In the simply laced case these
are `A_1`, `D_even`, `E_7`, `E_8`.
"""
function is_inc_type((X, n))
    X === :A && return n == 1
    X === :D && return iseven(n)
    X === :E && return n >= 7
    X === :I && return iseven(n)
    return X in (:B, :F, :G, :H)
end

# every component of the type t is an inc type and not one of bad_types
all_inc(t, bad_types) = all(x -> is_inc_type(x) && !(x in bad_types), t)

"""
    inc_subdiagrams(D, n; c_type = nothing, bad_types = (), subdiag = nothing)

The inc-type subdiagrams of rank `n` (spherical, and every component has an inc
type not in `bad_types`) as sorted vertex lists.  With `c_type`, exactly those of
that type instead; with `subdiag`, only those inside `subdiag`.
"""
function inc_subdiagrams(D::NamedTuple, n::Int; c_type = nothing, bad_types = (), subdiag = nothing)
    keep(t) = c_type === nothing ? all_inc(t, bad_types) : t == c_type
    return [I for (I, t) in spherical_subdiagrams(D, n; subdiag) if keep(t)]
end
inc_subdiagrams(labels::AbstractDict, n::Int; kwargs...) = inc_subdiagrams(diagram_data(labels), n; kwargs...)

"""
    inc_types_D(D, n; bad_types = (), subdiag = nothing)

The distinct types of the inc-type subdiagrams of rank `n`.
"""
function inc_types_D(D::NamedTuple, n::Int; bad_types = (), subdiag = nothing)
    return unique([t for (_, t) in spherical_subdiagrams(D, n; subdiag) if all_inc(t, bad_types)])
end
inc_types_D(labels::AbstractDict, n::Int; kwargs...) = inc_types_D(diagram_data(labels), n; kwargs...)

# ------------------------------------------------------------------------------
# 4. X_I
# ------------------------------------------------------------------------------

"""
    x_marker(D, I)

For an inc-type subdiagram `I`, `X_I` is the largest inc-type subdiagram of `I`
without `A_1` components, i.e. the union of the non-`A_1` components of `I`:
nothing can be chosen inside an `A_1` component, and any other inc component is
the largest choice inside itself.  These are exactly the vertices of `I` that
have a neighbour in `I`.
"""
x_marker(D, I::Vector{Int}) = [v for v in I if any(in(I), D.nbrs[v])]

# ------------------------------------------------------------------------------
# 5. The k-odd graph
# ------------------------------------------------------------------------------

"""
    k_odd(D, k; c_type = nothing, subdiag = nothing)

The k-odd graph `D^k` as `(g, verts, idx)`: vertex `i` of the graph `g` is the
subdiagram `verts[i]`, and `idx[verts[i]] == i`.  The vertices are the inc-type
subdiagrams of rank `k` (`c_type` and `subdiag` as in `inc_subdiagrams`), and
`I`, `J` are adjacent iff `X_I = X_J` and `J = (I \\ {s}) ∪ {t}` for an edge
`{s, t}` with odd `m_st`.

Instead of testing all pairs, the edges are constructed: since
`X_I = X_J ⊆ I ∩ J`, only a vertex `s ∈ I \\ X_I` can be replaced, by an odd
neighbour `t ∉ I`, and the resulting `J` is looked up in `idx`.
"""
function k_odd(D::NamedTuple, k::Int; c_type = nothing, subdiag = nothing)
    verts = inc_subdiagrams(D, k; c_type, subdiag)
    idx   = Dict(I => i for (i, I) in enumerate(verts))
    X     = [x_marker(D, I) for I in verts]

    edges = Vector{Int}[]
    for (i, I) in enumerate(verts), s in I
        s in X[i] && continue
        for t in D.oddnbrs[s]
            t in I && continue
            j = get(idx, sort!([filter(!=(s), I); t]), 0)
            j > i && X[j] == X[i] && push!(edges, [i, j])    # j > i: every edge once
        end
    end
    return graph_from_edges(Undirected, edges, length(verts)), verts, idx
end
k_odd(labels::AbstractDict, k::Int; kwargs...) = k_odd(diagram_data(labels), k; kwargs...)

"""
    odd_components(D, k; c_type = nothing, subdiag = nothing)

The connected components of `D^k` (see `k_odd`), each as a vector of subdiagrams.
"""
function odd_components(D::NamedTuple, k::Int; kwargs...)
    g, verts, _ = k_odd(D, k; kwargs...)
    return [verts[c] for c in connected_components(g)]
end
odd_components(labels::AbstractDict, k::Int; kwargs...) = odd_components(diagram_data(labels), k; kwargs...)
