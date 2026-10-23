using Random, StaticArrays


30 / 0.00001
a = round(2.9999999999999995e6)


# Vecteur standard
v = rand(3000000);

# NTuple de 1000 éléments (alternative immuable)
v_tuple = ntuple(_ -> rand(), 3000000);

# Vecteur de SVector{8} (1000 éléments)
v_svec = Vector{SVector{8, Float64}}(undef, 3000000)

v_svec2 = [SVector{7}(rand(7)) for _ in 1:a]

# Affichage des tailles en bits
println("Taille de v (Vector)         : ", Base.summarysize(v) * 8 / 1000000000, " Giga")
println("Taille de v_tuple (NTuple)   : ", Base.summarysize(v_tuple) * 8 / 1000000000, " Giga")
println("Taille de v_svec (Vector{SVector{8}}): ", Base.summarysize(v_svec) * 8 / 1000000000, " Giga")
println("Taille de v_svec2 (Vector{SVector{8}}): ", Base.summarysize(v_svec2) * 8 / 1000000000, " Giga")
