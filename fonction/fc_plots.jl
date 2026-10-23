#####################################################
#### Répertoire des fonction pour les graphiques ####
#####################################################

### Code générale ---- 
module module_graphique

# Fonction exporté 
export heatmap_compose

# Import packages
using CSV, Makie, CairoMakie, ColorSchemes, Plots, DataFrames  
using CairoMakie


### Création de la fonction de heatmap composé ----

# Paramètres 
#   nom_modele : "ParaLarve" ou "ParaNymphe" ou "ParaOeufs" ou "sansPara"
#               
#
#   axe1_name : string du nom du premier axe 
#               
#               
#   axe1_name: string du nom du premier axe 
#                     
#
#   nom_table : string du type de valeur affiché (moyenne, variance...)
#               


function heatmap_compose(iso, prov, nom_modele, DD, axe1_name, axe2_name,t_stud, mathematica = false)


    if Sys.KERNEL == :NT
        racine = "P:/"
    else
        racine = "/Users/Quentin/Desktop/"
    end
    
    # Charger les données de mathématica 
    if mathematica != false 
        data = CSV.read(""*racine*"Article1/tables/"*string(mathematica)*".csv", DataFrame; header=false)
    end 


    if prov == "MeSU"
        ## Charger des Dataframe
        resultat = CSV.read(""*racine*"Article1/tables/MeSU/"*nom_modele*"_heatmap_Densité_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultatStd = CSV.read(""*racine*"Article1/tables/MeSU/"*nom_modele*"_heatmap_Variance_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultaEq = CSV.read(""*racine*"Article1/tables/MeSU/"*nom_modele*"_Calque_Eq_type_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame; header=true)
        resultatFreq = CSV.read(""*racine*"Article1/tables/MeSU/"*nom_modele*"_heatmap_Frequence_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultatAir = CSV.read(""*racine*"Article1/tables/MeSU/"*nom_modele*"_heatmap_Integral_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
    
    elseif prov == "MeSU_low"
        ## Charger des Dataframe
        resultat = CSV.read(""*racine*"Article1/tables/MeSU_low/"*nom_modele*"_heatmap_Densité_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultatStd = CSV.read(""*racine*"Article1/tables/MeSU_low/"*nom_modele*"_heatmap_Variance_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultaEq = CSV.read(""*racine*"Article1/tables/MeSU_low/"*nom_modele*"_Calque_Eq_type_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame; header=true)
        resultatFreq = CSV.read(""*racine*"Article1/tables/MeSU_low/"*nom_modele*"_heatmap_Frequence_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultatAir = CSV.read(""*racine*"Article1/tables/MeSU_low/"*nom_modele*"_heatmap_Integral_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)

    elseif prov == "MeSUteste"
        ## Charger des Dataframe
        resultat = CSV.read(""*racine*"Article1/tables/MeSUteste/"*nom_modele*"_heatmap_Densité_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultatStd = CSV.read(""*racine*"Article1/tables/MeSUteste/"*nom_modele*"_heatmap_Variance_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultaEq = CSV.read(""*racine*"Article1/tables/MeSUteste/"*nom_modele*"_Calque_Eq_type_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame; header=true)
        resultatFreq = CSV.read(""*racine*"Article1/tables/MeSUteste/"*nom_modele*"_heatmap_Frequence_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultatAir = CSV.read(""*racine*"Article1/tables/MeSUteste/"*nom_modele*"_heatmap_Integral_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
    
    else 
        ## Charger des Dataframe
        resultat = CSV.read(""*racine*"Article1/tables/"*nom_modele*"/heatmap_Densité_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultatStd = CSV.read(""*racine*"Article1/tables/"*nom_modele*"/heatmap_Variance_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultaEq = CSV.read(""*racine*"Article1/tables/"*nom_modele*"/Calque_Eq_type_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame; header=true)
        resultatFreq = CSV.read(""*racine*"Article1/tables/"*nom_modele*"/heatmap_Frequence_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
        resultatAir = CSV.read(""*racine*"Article1/tables/"*nom_modele*"/heatmap_Integral_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", DataFrame)
    end

    # Extraire les axes 
    axe1 = names(resultat)[1:end-1]
    axe1 = tryparse.(Float64, axe1)
    axe2 = resultat[:,end]

    # Les passer en matrices 
    resultaEq = Matrix(resultaEq)
    resultatStd = Matrix(resultatStd)
    resultat = Matrix(resultat[:,1:end-1])
    resultatFreq = Matrix(resultatFreq)
    resultatAir = Matrix(resultatAir) 


    # Créer les masques
    masque0 = zeros(length(resultaEq[1,:]),length(resultaEq[:,1]))
    masque1 = zeros(length(resultaEq[1,:]),length(resultaEq[:,1]))
    masque2 = zeros(length(resultaEq[1,:]),length(resultaEq[:,1]))

    for i in 1:length(resultaEq[1,:])
        for j in 1:length(resultaEq[:,1])
            if resultaEq[i,j] == 0.0
                masque0[i,j] = 1.0
            else 
                masque0[i,j] = NaN
            end 

            if resultaEq[i,j] == 1.0
                masque1[i,j] = 1.0
            else 
                masque1[i,j] = NaN
            end 

            if resultaEq[i,j] == 2.0
                masque2[i,j] = 1.0
            else 
               masque2[i,j] = NaN
        end 
        end 
    end

   # Rescale les valeurs d'aire de la valeur pedant une période à la valeur par pas de temps 
   resultatAirRescale = replace!((resultatAir .* masque2) ./ (resultatFreq .* masque2), NaN => 0)
   
   # Définir les valeurs Min/Max
   Max = maximum(resultat)
   Min = minimum(resultat)

   MaxFreq = maximum(resultatFreq)
   MinFreq = minimum(resultatFreq)

   MaxAir = maximum(resultatAir)
   MinAir = Min

   MaxStd = maximum(resultatStd)
   MinStd = minimum(resultatStd)


   # Transposer la matrice pour contrer le biais d'affichage
   resultaEq = permutedims(resultaEq)
   resultat = permutedims(resultat)
   resultatStd = permutedims(resultatStd)
   resultatFreq = permutedims(resultatFreq)
   resultatAir = permutedims(resultatAir)
   masque0 = permutedims(masque0)
   masque1 = permutedims(masque1)
   masque2 = permutedims(masque2)


    # Création de la figure densité
    plt = Figure(resolution = (600, 600))   
    ax = Axis(plt[1,1], title="Heatmap compilé densité moyenne de larves")

    # Heatmap pour la seconde zone (Rouge)
    hm2 = CairoMakie.heatmap!(ax, axe1, axe2, resultat .* masque2, colormap=:Reds, colorrange=(Min, Max))

    # Heatmap pour la première zone (Bleu)
    hm1 = CairoMakie.heatmap!(ax, axe1, axe2, resultat .* masque1, colormap=:Blues, colorrange=(Min, Max))

    # Heatmap pour la troisième zone (Gris ou Bleu)
    hm0 = CairoMakie.heatmap!(ax, axe1, axe2, masque0, colormap=:Blues, colorrange=(Min, Max))

    # Ajout des isoclines (contours)
    if iso == true
        CairoMakie.contour!(ax, axe1, axe2, resultaEq, levels=2, color=:black, linewidth=1)
    end 

    # Fixation des limites des axes
limits!(ax, axe1[1], axe1[end], axe2[1], axe2[end])

    # Ajouter la limite de stabilité donnée par mathématica 
    if mathematica != false 
        x = data.Column1  
        y = data.Column2  
        lines!(ax, x, y, color=:black, linewidth=2)
    end 

    # Ajouter une Colorbar
    cb = Colorbar(plt[1,2], hm1, label="Densité moyenne")


   # Création de la figure fréquence
   plt2 = Figure(resolution = (600, 600))   
   ax = Axis(plt2[1,1], title="Heatmap compilé Période")

   # Heatmap pour la seconde zone
   hmFreq2 = CairoMakie.heatmap!(ax, axe1, axe2, resultatFreq .* masque2, colormap=:Reds, colorrange=(MinFreq,MaxFreq))
   hmFreq1 = CairoMakie.heatmap!(ax, axe1, axe2, resultatFreq .* masque1, colormap=:Blues, colorrange=(MinFreq,MaxFreq))
   hmFreq0 = CairoMakie.heatmap!(ax, axe1, axe2, resultatFreq .* masque0, colormap=:Blues, colorrange=(MinFreq,MaxFreq))
   # Ajouter une Colorbar
    cb = Colorbar(plt2[1,2], hmFreq2, label="durée de la Période")



   # Création de la figure Air sous la courbe
   pltAir = Figure(resolution = (600, 600))   
   ax = Axis(pltAir[1,1], title="Heatmap compilé Période")

   # Heatmap pour la seconde zone
   hmAir2 = CairoMakie.heatmap!(ax, axe1, axe2, resultatAir .* masque2, colormap=:Reds, colorrange=(MinAir,MaxAir))
   hmAir1 = CairoMakie.heatmap!(ax, axe1, axe2, resultat .* masque1, colormap=:Blues, colorrange=(MinAir,MaxAir))
   hmAir0 = CairoMakie.heatmap!(ax, axe1, axe2, resultat .* masque0, colormap=:Blues, colorrange=(MinAir,MaxAir))
   # Ajouter une Colorbar
    cb = Colorbar(pltAir[1,2], hmAir2, label="Air sous la courbe")



    # Création de la figure Écart type 
    pltStd = Figure(resolution = (600, 600))   
    ax = Axis(pltStd[1,1], title="Heatmap compilé écart-type des larves")

    # Heatmap pour la seconde zone (Rouge)
    hmStd = CairoMakie.heatmap!(ax, axe1, axe2, resultatStd, colormap=:Reds, colorrange=(MinStd, MaxStd))

    # Ajout des isoclines (contours)
    if iso == true
        CairoMakie.contour!(ax, axe1, axe2, resultaEq, levels=2, color=:black, linewidth=1)
    end 

    # Fixation des limites des axes
    limits!(ax, axe1[1], axe1[end], axe2[1], axe2[end])

    # Ajouter la limite de stabilité donnée par mathématica 
    if mathematica != false 
        x = data.Column1  
        y = data.Column2  
        lines!(ax, x, y, color=:black, linewidth=2)
    end 

    # Ajouter une Colorbar
    cb = Colorbar(pltStd[1,2], hmStd, label="Densité moyenne")


    save(""*racine*"Article1/graphs/"*nom_modele*"/heatmap_Combiné_Densité_"*DD*"_"*axe1_name*"_"*axe2_name*".png", plt)  
    save(""*racine*"Article1/graphs/"*nom_modele*"/heatmap_Combiné_Écart_type_"*DD*"_"*axe1_name*"_"*axe2_name*".png", pltStd)    
    save(""*racine*"Article1/graphs/"*nom_modele*"/heatmap_Combiné_Air_"*DD*"_"*axe1_name*"_"*axe2_name*".png", pltAir)
    save(""*racine*"Article1/graphs/"*nom_modele*"/heatmap_Combiné_Frequences_"*DD*"_"*axe1_name*"_"*axe2_name*".png", plt2)  
    save(""*racine*"Article1/graphs/temporaire/temp.png",pltStd) 
  

end












#### Création de la fonction de plot bifurcation + equlibre théorique

# Paramètres
#   chm : chemin d'acces au csv
#
#   bif : nom donnée bifurcation
#
#   dfp : (début, fin, pas)

function bifurcationReelle(chmMath, chmEmp, dfp)

    dataMath = CSV.read("/Users/Quentin/Desktop/Article1/tables/"*chmMath*".csv", DataFrame; header=false)
    dataEmps = CSV.read("/Users/Quentin/Desktop/Article1/tables/"*chmEmp*".csv", DataFrame; header=false)
    x = dataMath.Column1  
    y = dataMath.Column2  
    plt = plot(x, y, label="Données", xlabel="Colonne 1", ylabel="Colonne 2", title="Bifurcation et valeur d'équilibre")
    axe_range = range(dfp[1],dfp[2],dfp[3])
    for (i, axe) in enumerate(axe_range)
        scatter!(plt,
            fill(axe, length(dataEmps[i])),  # Toutes les valeurs de q (fixe) pour les résultats correspondants
            dataEmps[i],                    # Les premières valeurs extraites
            label=false,                               # Pas de légende
            color=:black,                              # Couleur noire
            alpha=0.4,                                 # Transparence faible (valeurs possibles : 0 à 1)
            ms=2)                                      # Taille des points (marker size)
    end
    display(plt)
    return dataMath
end

end
