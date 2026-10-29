# Sample computations.  From the repository folder, run
#
#     julia --project=. examples.jl
#
# or start `julia --project=.` and type `include("examples.jl")`.

include("isolated_subdiagrams.jl")

# --- 1. Isolated subdiagrams for diagrams from Diagrams.jl ---------------------

println("HE8:")
display(isolated_subdiagrams(labels_HE8, 9, affines_HE8, affines_par_HE8))
println("\nA1HE8:")
display(isolated_subdiagrams(labels_A1HE8, 10, affines_A1HE8, affines_par_A1HE8))

# The A1HE8E8 diagram has 24 vertices, so it is the slowest example; the first five
# ranks are quick, all 17 take longer.
println("\nA1HE8E8, ranks 1 to 5:")
display(isolated_subdiagrams(labels_A1HE8E8, 5, affines_A1HE8E8, affines_par_A1HE8E8))
# display(isolated_subdiagrams(labels_A1HE8E8, 17, affines_A1HE8E8, affines_par_A1HE8E8))

# --- 2. Types of subdiagrams ----------------------------------------------------

println("\nThe maximal affine subdiagrams of A1HE8E8:")
M = coxeter_matrix_from_labels(labels_A1HE8E8)
for (A, P) in zip(affines_A1HE8E8, affines_par_A1HE8E8)
    # coxeter_type returns `nothing` for infinite subgroups, e.g. affine ones
    println("affine subdiagram infinite: ", coxeter_type(M, A) === nothing,
            ",  spherical part: ", type_string(coxeter_type(M, P)))
end

# --- 3. Inc subdiagrams and k-odd graphs ----------------------------------------

println("\nInc types and k-odd graphs of HE8:")
D = diagram_data(labels_HE8)          # build once and reuse: it caches results
for k in 1:4
    println("rank $k inc types: ", join(type_string.(inc_types_D(D, k)), ", "))
end

t = [(:A, 1), (:A, 1)]                # type A₁²
g, verts, idx = k_odd(D, 2; c_type = t)
println("k-odd graph for ", type_string(t), ": ", length(verts), " vertices, ",
        length(odd_components(D, 2; c_type = t)), " connected component(s)")

