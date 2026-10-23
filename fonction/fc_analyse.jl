#######################################
#### Répertoir des Solveurs Maison ####
#######################################

### Code générale ---- 
module analyse

# Fonction exporté 
export raster2limite

# Import
using Plots, CSV, DataFrames

end



### Création de la fonction Solveur RK2 ----

# Paramètres 
#   matrix_path   :  chemin vers la matrice 
#
#   retard_vect : vecteur des retards
#
#   tspan : tuple de 2 Float64: 
#       tps au départ de simul
#       tps de fin de simul
#
#   u0  :  ajout initiale dans les pops 
#
#   pdt : Float64 du pas de temps utilisé
function raster2limite(matrix_path)

end


using Plots, CSV, DataFrames


matrix_path = "P:/Article1/tables/sansPara/heatmap_Variance_DD_q_deltaA.csv"

matrix_var = CSV.read(matrix_path, DataFrame)
#matrix_var = matrix_var[end:-1:1, :]

vect_abs = vec(names(matrix_var))[1:30]
vect_abs = round.(parse.(Float64, vect_abs),digits=1)
vect_ordo = vec(matrix_var[:,31])
matrix_var_0 = Matrix(matrix_var[:,1:30])

matrix_var_0 .= ifelse.(matrix_var_0 .< 1.1, 0, matrix_var_0)

matrix_lin = zeros(30,30)
lengthX = (length(matrix_var_0[1,:]))
lengthY = (length(matrix_var_0[:,1]))


for i in 1:lengthX
    for j in 1:lengthY
        if matrix_var_0[i,j] != 0  # La case courante n'est pas nulle
            # Vérifier les voisins
            if (i > 1 && matrix_var_0[i-1, j] == 0) ||  # Voisin en haut
               (i < lengthX && matrix_var_0[i+1, j] == 0) ||  # Voisin en bas
               (j > 1 && matrix_var_0[i, j-1] == 0) ||  # Voisin à gauche
               (j < lengthY && matrix_var_0[i, j+1] == 0)    # Voisin à droite
                matrix_lin[i,j] = matrix_var_0[i,j]  # Garder la valeur dans matrix_lin
            end
        end
    end
end


matrix_lin
matrix_lin_binaire = matrix_lin
matrix_lin_binaire .= ifelse.(matrix_lin_binaire .> 0, 1, matrix_lin_binaire)

p1 = heatmap( matrix_lin, color=:viridis)
p1 = heatmap(vect_abs, vect_ordo, matrix_lin, color=:viridis, title="Carte de stabilité")


















# Fonction pour trouver les lignes
# Fonction pour trouver les lignes en tenant compte des voisins diagonaux
function trouver_lignes(matrix)
    rows, cols = size(matrix)
    lignes = [] 

    # Fonction pour vérifier si un point est dans la matrice et égal à 1
    function est_valide(i, j)
        return i >= 1 && i <= rows && j >= 1 && j <= cols && matrix[i,j] == 1
    end

    # Liste des directions pour les voisins (haut, bas, gauche, droite, diagonales)
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
    # Matrice pour marquer les points déjà visités
    visité = falses(rows, cols)
    
    # Fonction pour effectuer une recherche en profondeur (DFS) pour trouver les voisins connectés
    function dfs(i, j, ligne)
        push!(ligne, (i, j))
        visité[i, j] = true
        
        for (di, dj) in directions
            ni, nj = i + di, j + dj
            if est_valide(ni, nj) && !visité[ni, nj]
                dfs(ni, nj, ligne)
            end
        end
    end

    # Parcourir chaque case de la matrice
    for i in 1:rows
        for j in 1:cols
            if matrix[i, j] == 1 && !visité[i, j]
                ligne = []  # Liste des points dans la ligne courante
                dfs(i, j, ligne)
                push!(lignes, ligne)  # Ajouter cette ligne trouvée à la liste des lignes
            end
        end
    end

    return lignes
end

# Exemple de matrice (10x10) avec des 1 et des 0
matrix = matrix_lin_binaire

# Trouver les lignes dans la matrice
lignes = trouver_lignes(matrix)

# Initialisation du plot
p = heatmap(matrix, color=:viridis, axis=false, title="Heatmap avec lignes");

# Superposer les lignes sur la heatmap
for ligne in lignes
    x_vals = [point[2] for point in ligne]  # Colonnes (abscisses)
    y_vals = [point[1] for point in ligne]  # Lignes (ordonnées)
    
    # Tracer les lignes reliant les points adjacents
    plot!(p, x_vals, y_vals, color=:red, linewidth=2, label="")
end

# Affichage du plot
display(p)