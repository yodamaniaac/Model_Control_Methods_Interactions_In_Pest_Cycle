# Obtenir le chemin absolu actuel
current_path = pwd()

# Diviser le chemin en parties
parts = splitdir(current_path)

# Recherche l'index 
index = findfirst(x -> occursin("Article1", x), parts)

# Couper après notre racine
if index !== nothing
    truncated_path = joinpath(parts[1:index]...)
    println("Chemin coupé après /Article1 : ", truncated_path)
else
    println("Le répertoire /Article1 n'est pas trouvé dans le chemin.")
end