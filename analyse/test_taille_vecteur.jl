using Random, StaticArrays



a = Int(round(100 / 0.0001))

# Vecteur standard
v = rand(a)

# NTuple de 1000 éléments (alternative immuable)
v_tuple = ntuple(_ -> rand(), a);

# Vecteur de SVector{8} (1000 éléments)
v_svec = Vector{SVector{8, Float64}}(undef, a)

v_svec2 = [SVector{7}(rand(7)) for _ in 1:a]

# Affichage des tailles en bits
println("Taille de v (Vector)         : ", Base.summarysize(v) * 8 / 1000000000, " Giga")
println("Taille de v_tuple (NTuple)   : ", Base.summarysize(v_tuple) * 8 / 1000000000, " Giga")
println("Taille de v_svec (Vector{SVector{8}}): ", Base.summarysize(v_svec) * 8 / 1000000000, " Giga")
println("Taille de v_svec2 (Vector{SVector{8}}): ", Base.summarysize(v_svec2) * 8 / 1000000000, " Giga")
