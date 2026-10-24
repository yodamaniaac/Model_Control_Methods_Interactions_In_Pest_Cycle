######################################################################################################
#### Répertoire des fonction pour la visualisation et l'analyse des fonction des modèles à retard ####
######################################################################################################

### Code générale ---- 
module module_analyse

# Fonction exporté 
export runDynamics, traceDynamics, heatMapStab, simulate_bifurcation, 
simulate_bifurcation_Final, plot_bifurcation_Final, Graph_survie_para_1D, 
without_parasitoid_model, eggs_parasitoid_model, larval_parasitoid_model, nymphal_parasitoid_model


# Import packages
using DifferentialEquations, Plots, NumericalIntegration, CSV, Trapz, FFTW, DelayDiffEq
using Trapz, Statistics, MAT, DelimitedFiles, Dates, MAT, DataFrames, ColorSchemes, OrdinaryDiffEq

# Final figures (as opposed to the scratch "graphs/temporaire" preview) are
# shared with the R comparison-table scripts under Supplementary_material/graphs.
const FINAL_GRAPH_DIR = joinpath(@__DIR__, "..", "graphs", "Final_Graph")
const BIFURCATION_GRAPH_DIR = joinpath(FINAL_GRAPH_DIR, "Bifurcation")
const FIGURE1D_GRAPH_DIR = joinpath(FINAL_GRAPH_DIR, "Figure1D")
mkpath(FINAL_GRAPH_DIR)
mkpath(BIFURCATION_GRAPH_DIR)
mkpath(FIGURE1D_GRAPH_DIR)

# Scratch preview used by printTemp()/the various plotting functions while a
# run is in progress -- always needs to exist locally, regardless of git state.
mkpath(joinpath(@__DIR__, "graphs", "temporaire"))
mkpath(joinpath(@__DIR__, "graphs", "sansPara"))


### Création de la fonction de dymanique ----

# Paramètres 
#   modele : fonction du modèle à retard 
#           
#   para   : tuple d'un tuple et d'un vecteur : 
#               tuple des paramètres du modèle 
#               vecteur des lags utilisé dans le modèle 
#
#   para_run : tuple contenant : 
#               vecteur des valeurs à 0
#               vecteur des valeurs historiques
#               tuple du début et fin de la simul
#               la méthode de résolution 
#               taile du pas de temps (facultatif)
#               
#   para_save: tuple des élément de localisation 
#               numéro du modèle 
#               Bool présence de densité dépendence sur les larves 
#               Bool save ou non 

        
function runDynamics(modele, para, para_run, para_save)
    
    # Extraire les paramètres de sauvgarde
    presence_dd = para_save[2]
    YN_save = para_save[3]
    nb_modele = para_save[4]

    # Extraire les paramètres du modèle 
    fixed_params = para[1]
    lags_fct =  para[2]
    ordre_para = para[3]

    p = fixed_params[ordre_para[1]]
    for (h) in 2:length(ordre_para)
        p = (p...,fixed_params[ordre_para[h]])
    end
 

    # Extraire les paramètres de simulation

    u0 = para_run[1](fixed_params,"u")
    h0 = para_run[1](fixed_params,"h")
    h(p,t) = h0
    tspan = para_run[2]
    alg = para_run[3]
    t_plot = para_run[4]
    if length(para_run) >= 5
        pdt = para_run[5]
    else 
        pdt = 0.0001
    end


    # Fixer l'horraire
    heure = string(now())

    # Lag total 
    lags = lags_fct(fixed_params)

    # Résoudre le problème selon le solveur choisis
    if alg == "fixe"
        # Définir le problème et le résoudre
        prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
        freq_save = 0.1

        sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=freq_save,abstol=1e-9, reltol=1e-9, maxiters = 1e9)

        # Extraire les solutions
        t_values = sol.t # Les temps où la solution est évaluée
        u_values = sol.u # Les valeurs des solutions correspondantes
    end

    t_values = collect(t_values)

    # Séléctionner les valeurs tracé
    if t_values[end] > t_plot
        t_last_values = [x for x in t_values if x > (t_values[end] - t_plot)]
        u_last_values = u_values[end-length(t_last_values)+1:end] 
    else 
        t_last_values = t_values
        u_last_values = u_values
    end 
    
    # Ne prendre qu'un point par unité de tps
    #t_lite = t_last_values[1:Int(round((1/pdt))):end]    Plus utile depuis l'ajout de REZ
    #u_lite = u_last_values[1:Int(round((1/pdt))):end]
    t_lite = t_last_values
    u_lite = u_last_values


    if YN_save == true
        
        if presence_dd == true
            DD = "DD"
        elseif presence_dd == false
            DD = "ssDD"
        end 

        nom_modele = para_save[1]
            writedlm( joinpath(@__DIR__, "tables")*"/"*nom_modele*"/dynamique_"*DD*"_"*heure*".csv", vcat(t_last_values', u_last_values'), ',')
    end
    return [u_lite, t_lite]
end


### Création de la fonction de traceDynamics ----

# Paramètres 
#   rez   : Vector{Vector{Any}} contenant les résultats des dynamiques 
#               
#   para_plot : tuple des paramètres pour l'affichage des courbes
#               tuple du numéro des variables tracé
#               tuple du nom des varaibles
#               taille du segment de la dynamique tracé

function traceDynamics(rez, para_plot)

    plot_title = para_plot[1]
    x_name = para_plot[2]
    y_name = para_plot[3]
    curve_plot = para_plot[4]
    labels = para_plot[5]

    # Créer les couleurs 
    palette = ColorSchemes.colorschemes[:twelvebitrainbow]
    positions = collect(range(0, 1; length=12))
    indices = [4, 6, 10, 2]
    colors = [RGB{Float64}(get(palette, positions[i])) for i in indices]

    # Créer le plot vide
    plt = plot(title=plot_title, xlabel = x_name, ylabel = y_name, lw=4, grid=false, titlefont = font(15), legend=:topright, legendfont = font(10),yscale =:log10)

    # Ajouter les courbes
    for i in 1:length(rez)
            plot!(plt, rez[i][2],  [u[curve_plot] + 1 for u in rez[i][1]], 
            label=labels[i], color=colors[i], linewidth = 2)
    end

    display(plt)  
    
    return(plt)
end 


### Création de la fonction heatMapStab ----

# Paramètres 
#   para   : tuple d'un tuple et d'un vecteur : 
#               dictionnaire des paramètres du modèle et de leur valeur
#               vecteur des lags utilisé dans le modèle 
#               ordre des para dans p
#
#   para_run : tuple contenant : 
#               fonction créant les valeurs des varialbes à 0 et en historique
#               tuple du début et fin de la simul
#               la méthode de résolution 
#      
#   para_matrice : tuple contenant : 
#               tuple des deux noms de paramètres qui serviront d'axe à la matrice
#               tuple des valeur min et max des deux paramètres               
#               tuple des dimensions des deux axes de la matrice
#               int des X dernier pas de temps étudié
#               
#         
#   para_save: tuple des élément de localisation 
#               numéro du modèle 
#               Bool présence de densité dépendence sur les larves 
#               Bool save ou non 
#               
#
#   para_plot : tuple des paramètres pour l'affichage des courbes
#               tuple du numéro de la variables dont la stabilité est observé
#               tracer la heatmap de stabilité
#               tracer la heatmap des moyennes 

function heatMapStab(modele,para,para_run,para_matrice,para_save,para_plot)

    # Extraire les paramètres du modèle 
    fixed_params = para[1]
    lags_fct =  para[2]
    ordre_para = para[3]

    # Extraire les paramètres de simulation
    tspan = para_run[2]
    alg = para_run[3]
    pdt = para_run[4]
    

    # Etraire les paramètre pour la création de la heapmap
    axe1_name = string(para_matrice[1][1])
    axe2_name = string(para_matrice[1][2])
    axe1_range = range(para_matrice[2][1][1],para_matrice[2][1][2],para_matrice[3][1])
    axe2_range = range(para_matrice[2][2][1],para_matrice[2][2][2],para_matrice[3][2])

    resultsVar = zeros(para_matrice[3][1], para_matrice[3][2])
    resultMean = zeros(para_matrice[3][1], para_matrice[3][2])
    resultsMax = zeros(para_matrice[3][1], para_matrice[3][2])
    resultsMin = zeros(para_matrice[3][1], para_matrice[3][2])
    resultsFreq = zeros(para_matrice[3][1], para_matrice[3][2])
    resultsIntegral = zeros(para_matrice[3][1], para_matrice[3][2])

    t_stud = para_matrice[4]

    # Extraire les paramètres de sauvgarde
    nb_modele = para_save[1]
    presence_dd = para_save[2]
    YN_save = para_save[3]

    # Extraire les paramètres pour le graphique
    pop_obs = para_plot[1]
    YN_stab = para_plot[2]
    YN_mean = para_plot[3]

    # Boucle créatrice de la heatmap
    for (i, axe1) in enumerate(axe1_range)
        for (j, axe2) in enumerate(axe2_range)
           
            # Créer les paramètres dynamiques pour ce run
            p = Dict(axe1_name => axe1, axe2_name => axe2)
            for (key, val) in fixed_params
                p[key] = val
            end

            # Transformer les dico en tuple
            lags = lags_fct(p)

            p_tuple = p[ordre_para[1]]
            for (h) in 2:length(ordre_para)
                p_tuple = (p_tuple...,p[ordre_para[h]])
            end
            p_dico = p
            p = p_tuple 


            # Définir les valeur des varaibles en 0 et leur historique
            u0 = para_run[1](p_dico,"u")
            h0 = para_run[1](p_dico,"h")
            h(p,t) = h0

            # Lag total 
            totLag = p_dico["tauL"] + p_dico["tauL"] + p_dico["tauL"]

            # Résoudre le problème selon le solveur choisis
            if alg == "fixe"
                # Définir le problème 
                frq_save = 0.001
                prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)

                # Définir une condition d'arrêt basée sur l'écart type
                function condition(u, t, integrator)
                    # Vérifier toutes les 200 unités de temps
                    if t >= 200.0*2.5 && t % 200.0 < integrator.dt
                        # Calculer l'écart type des valeurs de u durant les 200 dernières unités de temps

                        # Filtrez les valeurs qui sont supérieures ou égales au seuil
                        recent_time = integrator.sol.t[integrator.sol.t .>= integrator.sol.t[end]-200]
                        recent_values = [v[1] for v in integrator.sol.u[end-(length(recent_time)-1):end]]
                        std_u = std(recent_values)

                        if std_u < 0.01
                            return true
                        elseif (abs(maximum(recent_values[Int(round(length(recent_values)/2)):end]) - 
                            maximum(recent_values[1:Int(round(length(recent_values)/2))-1]))) / mean(recent_values) < 0.01 &&
                            abs(minimum(recent_values[Int(round(length(recent_values)/2)):end]) - 
                            minimum(recent_values[1:Int(round(length(recent_values)/2))-1])) / mean(recent_values) < 0.01
                            return true 
                        end 
                    end
                    return false
                end

                # Affecter la condition d'arrêt
                function affect!(integrator)
                    terminate!(integrator)
                end

                # Définir un callback pour arrêter l'intégration
                cb = DiscreteCallback(condition, affect!, save_positions=(false, false))

                sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=frq_save,abstol=1e-9, reltol=1e-9, callback=cb) 
                #sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=frq_save) 

       
                # Extraire les solutions
                t_values = sol.t # Les temps où la solution est évaluée
                u_values = sol.u # Les valeurs des solutions correspondantes
                eq_type = 1
            else 
                # Définir le problème et le résoudre
                prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
                sol = solve(prob, MethodOfSteps(Tsit5()))  
     
                # Extraire les solutions
                t_values = sol.t # Les temps où la solution est évaluée
                u_values = sol.u # Les valeurs des solutions correspondantes
            end
            
            #println("")
            #println("para p="*string(p_tuple[1])*" deltaU="*string(p_tuple[4])*"")
            #println("dynamique faite : t="*string(t_values[end])*" u="*string(u_values[end])*" eq_type="*string(eq_type)*"")
            
            # Sélectionner les dernières unités de temps 
            t_values = collect(t_values)

            if t_values[end] > t_stud
                t_last_values = [x for x in t_values if x > (t_values[end] - t_stud)]
                #u_last_values = u_values[end-(t_stud-1):end] 
                u_last_values = u_values[end-length(t_last_values)+1:end] 

            else 
                t_last_values = t_values
                u_last_values = u_values
            end


            # Selectionner la population observé
            u_obs = [v[pop_obs] for v in u_last_values]

            # Calculer la variance et la moyenne des 100 derniers pas de temps
            resultMean[j, i] = mean(u_obs)
            resultsVar[j, i] = std(u_obs)
            resultsMax[j, i] = maximum(u_obs)
            resultsMin[j, i] = minimum(u_obs)

        end
    end

    # Nettoyage des données invalides
    resultsVar[resultMean .== 0] .= 0

    # Tracer les heatmaps
    if YN_stab == true & YN_mean == false
        plt = heatmap(axe1_range, axe2_range, resultsVar, color=:viridis,
            xlabel=axe1_name, ylabel=axe2_name, title="Carte de stabilité")
            type_heatmap = "stabilité"

    elseif  YN_stab == false & YN_mean == true
        plt = heatmap(axe1_range, axe2_range, resultMean, color=:viridis,
            xlabel=axe1_name, ylabel=axe2_name, title="Carte des valeurs à l'Équilibre")
            type_heatmap = "moyenne"

    elseif YN_stab == true & YN_mean == true
        p1 = heatmap(axe1_range, axe2_range, resultsVar, color=:viridis,
            xlabel=axe1_name, ylabel=axe2_name, title="Carte de stabilité")

        p2 = heatmap(axe1_range, axe2_range, resultMean, color=:viridis,
            xlabel=axe1_name, ylabel=axe2_name, title="Carte des valeurs à l'Équilibre")

        plt = plot(p1, p2, layout=(1, 2), size=(1000, 400))
        type_heatmap = "complète"

    end
    #display(plt)

    if YN_save == true
         
        # Fixer les paras de sauvgarde
        heure = string(now())
        nom_modele = ["sansPara","ParaOeufs","ParaLarve","ParaNymphe"][nb_modele+1]

        if presence_dd == true
            DD = "DD"
        elseif presence_dd == false
            DD = "ssDD"
        end 
        
        # Créer le Dataframe des données à sauvgarder 
        
        # Convertir la matrice en DataFrame
        resultsVar = DataFrame(resultsVar, Symbol.(axe1_range))
        resultMean = DataFrame(resultMean, Symbol.(axe1_range))
        resultsMax = DataFrame(resultsMax, Symbol.(axe1_range))
        resultsMin = DataFrame(resultsMin, Symbol.(axe1_range))        
        resultsFreq = DataFrame(resultsFreq, Symbol.(axe1_range))
        resultsIntegral = DataFrame(resultsIntegral, Symbol.(axe1_range))

        # Ajouter une colonne pour les légendes des lignes
        resultsVar[!, "Row Labels"] = axe2_range
        resultMean[!, "Row Labels"] = axe2_range


        nom_modele = ["sansPara","ParaOeufs","ParaLarve","ParaNymphe"][nb_modele+1]
        if Sys.KERNEL == :NT
            savefig(plt, "$(@__DIR__)/graphs/temporaire/temp.png") 
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Variance_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsVar) 
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Densité_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultMean)
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Maximum_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsMax;)
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Minimum_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsMin;)
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Frequence_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsFreq;)
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Integral_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsIntegral;)
            savefig(plt, joinpath(FINAL_GRAPH_DIR, nom_modele*"_heatmap_"*type_heatmap*"_"*DD*"_"*axe1_name*"_"*axe2_name*".png"))      
        else
            savefig(plt, "$(@__DIR__)/graphs/temporaire/temp.png") 
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Variance_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsVar) 
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Densité_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultMean)
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Maximum_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsMax;)
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Minimum_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsMin;)
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Frequence_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsFreq;)
            CSV.write("$(@__DIR__)/tables/"*nom_modele*"/heatmap_Integral_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsIntegral;)
            savefig(plt, joinpath(FINAL_GRAPH_DIR, nom_modele*"_heatmap_"*type_heatmap*"_"*DD*"_"*axe1_name*"_"*axe2_name*".png"))    
        end
    elseif YN_save == false 
        if Sys.KERNEL == :NT
            savefig(plt, "$(@__DIR__)/graphs/sansPara/temp.png") 
        else
            savefig(plt, "$(@__DIR__)/graphs/sansPara/temp.png") 
        end    
    end
    return resultMean,resultsVar
end



### Création de la fonction diagramme de bifurcation ----

# Paramètres 
#   para   : tuple d'un tuple et d'un vecteur : 
#               dictionnaire des paramètres du modèle et de leur valeur
#               vecteur des string des lags utilisé dans le modèle 
#               ordre des para dans p
#
#   para_run : tuple contenant : 
#               fonction créant les valeurs des varialbes à 0 et en historique
#               tuple du début et fin de la simul
#               la méthode de résolution 
#               taile du pas de temps (facultatif)
#      
#   para_bifucation : tuple contenant : 
#               noms du paramètres qui servira d'axeà la bifurcation
#               min et max du paramètre              
#               nb d'itération de calcule 
#               int des X dernier pas de temps étudié
#          
#         
#   para_save : tuple des élément de localisation 
#               numéro du modèle 
#               Bool présence de densité dépendence sur les larves 
#               Bool save ou non 
#               
#
#   plot : tuple des paramètres pour l'affichage des courbes
#               tuple du numéro de la variables dont la stabilité est observé
#               affiche toutes les valeurs prise ou juste les max/min 
#
#   solveur_type : String type de solveur ("default", "fixe", "fixeSimple")
function simulate_bifurcation(modele, para, para_run, para_bifucation, para_save, para_plot)
    
    # Extraire les paramètres du modèle 
    fixed_params = para[1]
    lags_fct =  para[2]
    ordre_para = para[3]

    # Extraire les paramètres de simulation
    tspan = para_run[2]
    alg = para_run[3]
    if length(para_run) >= 4
        pdt = para_run[4]
    else 
        pdt = 0.01
    end
    

    # Etraire les paramètre pour la création de la bifurcation
    axe_name = string(para_bifucation[1])
    axe_range = range(para_bifucation[2][1],para_bifucation[2][2],para_bifucation[3])
    t_stud = para_bifucation[4]

    # Préparation du vecteur contenant les résultats 
    bifurcation_rez = []

    # Extraire les paramètres de sauvgarde
    nb_modele = para_save[1]
    presence_dd = para_save[2]
    YN_save = para_save[3]

    # Extraire les paramètres pour le graphique
    pop_obs = para_plot[1]
    plot_all = para_plot[2]


    # Boucle créatrice de la bifurcation
    for (i, axe) in enumerate(axe_range)

        
        # Créer les paramètres dynamiques pour ce run
        p = Dict()
        
        for (key, val) in fixed_params
            p[key] = val
        end

        p[axe_name] = axe

        # Transformer les dico en tuple
        lags = lags_fct(p)

        p_tuple = p[ordre_para[1]]
        for (h) in 2:length(ordre_para)
            p_tuple = (p_tuple...,p[ordre_para[h]])
        end
        p_dico = p
        p = p_tuple 


        # Définir les valeur des varaibles en 0 et leur historique
        u0 = para_run[1](p_dico,"u")
        h0 = para_run[1](p_dico,"h")
        h(p,t) = h0

        # Lag total 
        totLag = p_dico["tauL"] + p_dico["tauL"] + p_dico["tauL"]

        # Résoudre le problème selon le solveur choisis
        if alg == "fixe"
              # Définir le problème et le résoudre
              prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
              if pdt <= 0.01
                freq_save = 0.01
              else
                freq_save = pdt
              end
              #sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=freq_save)  
              sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=freq_save,abstol=1e-9, reltol=1e-9, maxiters = 1e9)
       
              # Extraire les solutions
              t_values = sol.t # Les temps où la solution est évaluée
              u_values = sol.u # Les valeurs des solutions correspondantes
        else 
            # Définir le problème et le résoudre
            prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
            sol = solve(prob, MethodOfSteps(Tsit5()))  
     
            # Extraire les solutions
            t_values = sol.t # Les temps où la solution est évaluée
            u_values = sol.u # Les valeurs des solutions correspondantes
        end

        # Sélectionner les dernières unités de temps 
        t_values = collect(t_values)
            if t_values[end] > t_stud
            t_last_values = [x for x in t_values if x > (t_values[end] - t_stud)]
            u_last_values = u_values[end-length(t_last_values)+1:end] 
        else 
            t_last_values = t_values
            u_last_values = u_values
        end
        
        # Ne prendre qu'un point par unité de tps
        t_lite = t_last_values[1:Int(round((1/pdt))):end]
        u_lite = u_last_values[1:Int(round((1/pdt))):end]

        # Afficher tout ou min/max
        if plot_all == true
            push!(bifurcation_rez, [u_lite])
        elseif nb_modele == 0
            MaxMinMean = [[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
            for ops in 1:4
                tot_ops = [v[ops] for v in u_lite] 
                MaxMinMean[1][ops] = maximum(tot_ops)
                MaxMinMean[2][ops] = minimum(tot_ops)
                MaxMinMean[3][ops] = mean(tot_ops)
            end 
            push!(bifurcation_rez, MaxMinMean)
        else 
            MaxMinMean = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]
            for ops in 1:7
                tot_ops = [v[ops] for v in u_lite] 
                MaxMinMean[1][ops] = maximum(tot_ops)
                MaxMinMean[2][ops] = minimum(tot_ops)
                MaxMinMean[3][ops] = mean(tot_ops)
            end 
            push!(bifurcation_rez, MaxMinMean)
        end 


    end


    # Tracer le diagramme de bifurcation
    

    if plot_all == true
        yname="Valeurs asymptotiques"
        drk = 0.04

        plt = plot(title="Diagramme de bifurcation Modèle "*string(nb_modele)*" ",
        xlabel=axe_name,
        ylabel=yname)

        # On réduit bifurcation_rezz
        reduced_bifurcation = [[[[v[pop_obs]] for v in subsub] for subsub in sub] for sub in bifurcation_rez]

        for (i, axe) in enumerate(axe_range)
            scatter!(plt,
                fill(axe, length(reduced_bifurcation[i])),  # Toutes les valeurs de q (fixe) pour les résultats correspondants
                reduced_bifurcation[i],                    # Les premières valeurs extraites
                label=false,                               # Pas de légende
                color=:black,                              # Couleur noire
                alpha=drk,                                 # Transparence faible (valeurs possibles : 0 à 1)
                ms=2)                                      # Taille des points (marker size)
        end

    elseif nb_modele == 0
        yname="Minimum et Maximim à l'asymptote"
        drk = 1

        plt = plot(title="Diagramme de bifurcation Modèle "*string(nb_modele)*" ",
        xlabel=axe_name,
        ylabel=yname)

        colors = [:green, :blue, :purple, :brown]
        labels = ["Oeufs","Larves","Nymphes","Adultes"]

        reduced_bifurcation = bifurcation_rez

        for (i, axe) in enumerate(axe_range)
            for ii in 1:2
                for iii in 2:length(bifurcation_rez[i][ii])-1
                    point = bifurcation_rez[i][ii][iii]
                    scatter!(plt, [axe], [point], label=false, color=colors[iii],
                            alpha=drk, ms=4)
                end
            end 
            
        end

        # Tracer les lignes
        vect_line = [v[3] for v in bifurcation_rez] 
        x_line = axe_range
        for i in 2:length(vect_line[1])-1
            y_line = [v[i] for v in vect_line] 
            plot!(plt, x_line, y_line, linewidth=1, color=colors[i], label = labels[i])
        end

        # Appliquer les limites de l'axe Y après avoir tout tracé
        #plot!(plt, ylims=(0, 3500))

    else
        yname="Minimum et Maximim à l'asymptote"
        drk = 1

        plt = plot(title="Diagramme de bifurcation Modèle "*string(nb_modele)*" ",
        xlabel=axe_name,
        ylabel=yname)

        colors = [:green, :blue, :purple, :brown, :orange, :red]
        labels = ["Oeufs","Larves","Nymphes","Adultes","Parasitoïdes juvéniles", "Parasitoïde adulte"]

        reduced_bifurcation = bifurcation_rez

        for (i, axe) in enumerate(axe_range)
            for ii in 1:2
                for iii in 3:length(bifurcation_rez[i][ii])-1
                    point = bifurcation_rez[i][ii][iii]
                    scatter!(plt, [axe], [point], label=false, color=colors[iii],
                            alpha=drk, ms=4)

                     plot!(
                        plt,
                        axe_range,
                        max_vals,
                        fillrange = min_vals,
                        label = false,
                        color = colors[o],
                        alpha = 0.3,
                        linewidth = 0,
                        fillalpha = 0.15,
                    )
                end
            end 
            
        end

        # Tracer les lignes
        vect_line = [v[3] for v in bifurcation_rez] 
        x_line = axe_range
        for i in 3:length(vect_line[1])-1
            y_line = [v[i] for v in vect_line] 
            plot!(plt, x_line, y_line, linewidth=1, color=colors[i], label = labels[i])
        end

        # Appliquer les limites de l'axe Y après avoir tout tracé
        plot!(plt, ylims=(0, 200))

    end


    #display(plt)

             
    #Fixer les paras de sauvgarde
    heure = string(now())
    nom_modele = ["sansPara","ParaOeufs","ParaLarve","ParaNymphe"][nb_modele+1]

    if presence_dd == true
        DD = "DD"
    elseif presence_dd == false
        DD = "ssDD"
    end 

    if plot_all == true
        pltAll = "all"
    else
        pltAll = "MinMax"
    end
        
    nom_modele = ["sansPara","ParaOeufs","ParaLarve","ParaNymphe"][nb_modele+1]

    if YN_save == true

        if Sys.KERNEL == :NT
            savefig(plt, "$(@__DIR__)/graphs/temporaire/temp.png") 
            #writedlm("$(@__DIR__)/tables/"*nom_modele*"/dataBifurcation_"*DD*"_"*pltAll*"_"*heure*".csv", bifurcation_rez, ',') 
            #savefig(plt, joinpath(FINAL_GRAPH_DIR, nom_modele*"_dataBifurcation_"*DD*"_"*pltAll*"_"*heure*".png"))      
        else
            writedlm("$(@__DIR__)/tables/"*nom_modele*"/dataBifurcation_"*DD*"_"*pltAll*"_"*heure*".csv", bifurcation_rez, ',') 
            savefig(plt, joinpath(FINAL_GRAPH_DIR, nom_modele*"_dataBifurcation_"*DD*"_"*pltAll*"_"*heure*".png")) 
            savefig(plt, "$(@__DIR__)/graphs/temporaire/temp.png")      
        end
    elseif YN_save == false 
        if Sys.KERNEL == :NT
            savefig(plt, "$(@__DIR__)/graphs/temporaire/temp.png") 
        else
            writedlm("$(@__DIR__)/tables/"*nom_modele*"/dataBifurcation_"*DD*"_"*pltAll*"_"*heure*".csv", bifurcation_rez, ',') 
            savefig(plt, "$(@__DIR__)/graphs/temporaire/temp.png") 
        end    
    end

    return reduced_bifurcation
end







### Création de la fonction Affichant le graphique temporaire ----

function printTemp(OnOff=true)
    if OnOff == true
        if Sys.KERNEL == :NT
            run(`cmd /c start "$(@__DIR__)/graphs/temporaire/temp.png"`)       
        else
            run(`open "$(@__DIR__)/graphs/temporaire/temp.png"`)
        end
    end 
end






### Création de la fonction diagramme de bifurcation Final  ----

# Paramètres 
#
#   modele : fonction : fonction permettant la simulation numérique des différents modèles
#
#   para   : tuple d'un tuple et d'un vecteur : 
#               dictionnaire des paramètres du modèle et de leur valeur
#               vecteur des string des lags utilisé dans le modèle 
#               ordre des para dans p
#
#   para_run : tuple contenant : 
#               fonction créant les valeurs des varialbes à 0 et en historique
#               tuple du début et fin de la simul
#      
#   para_bifucation : tuple contenant : 
#               noms des paramètres qui serviront d'axe à la bifurcation
#               min et max du paramètre              
#               nb d'itération de calcule 
#               int des X dernier pas de temps étudié
#          
#         
#   para_save : tuple des élément de localisation 
#               numéro du modèle 
#               Bool save ou non 
#               Save Y or N
#   max : valeur maximum prose par la simulation
#                  
#   type : string du type de bifurcation étudié
#
#   solveur_type : String type de solveur ("default", "fixe", "fixeSimple")
function simulate_bifurcation_Final(modele, para, para_run, para_bifucation, para_save, max, type, popobs=2)

    if max isa Float64 
    else
        max=1
    end

    # Extraire les paramètres du modèle 
    fixed_params = para[1]
    lags_fct =  para[2]
    ordre_para = para[3]

    # Extraire les paramètres de simulation
    tspan = para_run[2]
    if length(para_run) >= 3
        pdt = para_run[3]
    else 
        pdt = 0.0001
    end

    # Extraire les paramètres de sauvgarde
    nb_modele = para_save[4]
    YN_save = para_save[3]
    
    # Extraire les paramètres pour la création de la bifurcation
     axe_range = range(para_bifucation[2][1],para_bifucation[2][2],para_bifucation[3])


    t_stud = para_bifucation[4]


    # Préparation du vecteur contenant les résultats totaux
    bifurcation_rez_tot = Vector{Vector{Vector{Float64}}}()
    bifurcation_rez = Vector{Vector{Float64}}()

    #for (o,morta) in enumerate(para_bifucation[1])
    for (o,morta) in enumerate(para_bifucation[5])

        if type[1] == "1D"
            if morta == "Without Mortality"
                fixed_params["deltaU"]=0
                fixed_params["deltaL"]=0
                fixed_params["deltaM"]=0 
            elseif morta == "Egg Mortality"
                fixed_params["deltaU"]=0.1
                fixed_params["deltaL"]=0
                fixed_params["deltaM"]=0
            elseif morta == "Larvae Mortality"
                fixed_params["deltaU"]=0
                fixed_params["deltaL"]=0.1
                fixed_params["deltaM"]=0
            elseif morta == "Nymphae Mortality"
                fixed_params["deltaU"]=0
                fixed_params["deltaL"]=0
                fixed_params["deltaM"]=0.1
            end
        end 
        
        # Fixer le nom des axes 
        #axe_name = string(morta)
        if type[1] == "bifurcation_mortalite"
            axe_name = string(para_bifucation[1][o])
        else 
            axe_name = string(para_bifucation[1])
        end 

        # Préparation du vecteur contenant les résultats locaux
        bifurcation_rez = Vector{Vector{Float64}}()


        # Boucle créatrice de la bifurcation
        for (i, axe) in enumerate(axe_range)
            
            # Créer les paramètres dynamiques pour ce run
            p = Dict()
            
            for (key, val) in fixed_params
                p[key] = val
            end

            p[axe_name] = axe

            # Transformer les dico en tuple
            lags = lags_fct(p)

            p_tuple = p[ordre_para[1]]
            for (h) in 2:length(ordre_para)
                p_tuple = (p_tuple...,p[ordre_para[h]])
            end
            p_dico = p
            p = p_tuple 

            # Définir les valeur des varaibles en 0 et leur historique
            u0 = para_run[1](p_dico,"u")
            h0 = para_run[1](p_dico,"h")
            h(p,t) = h0


            # Définir le problème et le résoudre
            prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
            freq_save = 0.1

            sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=freq_save,abstol=1e-9, reltol=1e-9, maxiters = 1e9)
        
            # Extraire les solutions
            t_values = sol.t # Les temps où la solution est évaluée
            u_values = sol.u # Les valeurs des solutions correspondantes
            

            # Sélectionner les dernières unités de temps 
            t_values = collect(t_values)
                if t_values[end] > t_stud
                t_last_values = [x for x in t_values if x > (t_values[end] - t_stud)]
                u_last_values = u_values[end-length(t_last_values)+1:end] 
            else 
                t_last_values = t_values
                u_last_values = u_values
            end
            
            # Ne prendre qu'un point par unité de tps
            #u_lite = u_last_values[1:Int(round((1/pdt))):end]

            # Garder min/max/mean
            MaxMinMean = [0.0,0.0,0.0,0.0]
            tot_ops = [v[popobs] for v in u_last_values]
            MaxMinMean[1] = maximum(tot_ops)/max
            MaxMinMean[2] = minimum(tot_ops)/max
            MaxMinMean[3] = mean(tot_ops)/max

            # Ajouter la valeur moyenne du para 
            if nb_modele != 0 && mean([v[6] for v in u_last_values]) > 0.001
                MaxMinMean[4] = 1.0 
            end

            # Inserer les points pour cette mortalité dans la sortie total
            push!(bifurcation_rez, MaxMinMean)
                 
        end  
        push!(bifurcation_rez_tot, bifurcation_rez)
    end

    return bifurcation_rez_tot
end


### Création de la fonction diagramme de bifurcation Final  ----

# Paramètres 
#   bifurcation_rez_tot   : données utilisé pour le graph
#      
#   para_bifucation : tuple contenant : 
#               noms des paramètres qui serviront d'axe à la bifurcation
#               min et max du paramètre              
#               nb d'itération de calcule 
#               int des X dernier pas de temps étudié
#          
#         
#   para_save : tuple des élément de localisation 
#               numéro du modèle 
#               Bool save ou non 
#               Save Y or N

function plot_bifurcation_Final(bifurcation_rez_tot, para_bifucation, para_save, para_plot_panel, para_plot)

    # Extraire les paramètres de sauvgarde
    nom_modele = para_save[1]
    YN_save = para_save[3]
    nb_modele = para_save[4]

    # Extraire les paramètres pour la création de la bifurcation
    axe_range = range(para_bifucation[2][1],para_bifucation[2][2],para_bifucation[3])



    # Tracer le diagramme de bifurcation
    plt = plot(title="",
    xlabel=para_plot[2],
    ylabel=para_plot[3],
    grid=false)

 
    if para_bifucation[6] == "log"
        Plots.plot!(yticks = (log10.([1,10,100,1000,10000]),["1","10","100","1000","10000"]))
    elseif para_bifucation[6] == "ratio"
        Plots.plot!(yticks = ([0,20,40,60,80,100],["0%","20%","40%","60%","80%","100%"])) 
        end 

    if para_plot[1] isa Vector{Float64}
        Plots.plot!(ylims=(para_plot[1][1],para_plot[1][2]))
    end 

    labels = para_bifucation[5]
    palette = ColorSchemes.colorschemes[:twelvebitrainbow]
    positions = collect(range(0, 1; length=12))
    indices = [2, 4, 6, 10]
    colors = [RGB{Float64}(get(palette, positions[i])) for i in indices]


    #for (o,morta) in enumerate(para_bifucation[1])
    for (o,morta) in enumerate(para_bifucation[5])

        # Afficher la zone min/max,
        min_vals = [bifurcation_rez_tot[o][i][1] for (i, _) in enumerate(axe_range)]
        max_vals = [bifurcation_rez_tot[o][i][2] for (i, _) in enumerate(axe_range)]

        # Ajuster les valeur en log pour éviter les -inf
        if para_bifucation[6] == "log"
            min_vals = log10.(min_vals .+ 1)
            max_vals = log10.(max_vals .+ 1)
        elseif para_bifucation[6] == "ratio"
            min_vals = exp.(min_vals)*100
            max_vals = exp.(max_vals)*100
        end 


        Plots.plot!(
            plt,
            axe_range,
            max_vals,
            fillrange = min_vals,
            label = para_bifucation[5][o],
            color = colors[o],
            alpha = 0.6,
            linewidth = 0,
            fillalpha = 0.15,
        )

        if para_bifucation[6] == "ratio"
            Plots.plot!(
                legend = false
            )
        end 


        # Tracer les lignes
        vect_line = [v[3] for v in bifurcation_rez_tot[o]] 
        
        # Ajuster les valeur en log pour éviter les -inf
        if para_bifucation[6] == "log"
            vect_line = log10.(vect_line .+ 1) 
        elseif para_bifucation[6] == "ratio"
            vect_line = exp.(vect_line)*100 
        end 

        presence_para = [v[4] for v in bifurcation_rez_tot[o]] 
        if para_bifucation[6] == "ratio"
            presence_para = ones(length(presence_para))
        end 
        x_line = axe_range

        # Trouver l'indice de changement
        change_idx = findfirst(presence_para .!= presence_para[1])
        #if isa(findfirst(presence_para .!= presence_para[1]),Nothing) !=  true
        #    change_idx = change_idx = findfirst(presence_para .!= presence_para[1])-1
        #end

        if isnothing(change_idx)
            # Pas de changement
            style = presence_para[1] == 0.0 ? :dashdotdot : :solid
            plot!(x_line, vect_line, linestyle=style, linewidth=2, color=colors[o], label = false)
        else
            # Première partie, jusqu'au point juste avant le changement
            plot!(x_line[1:change_idx-1], vect_line[1:change_idx-1],
                linestyle=presence_para[1] == 0.0 ? :dashdotdot : :solid, linewidth=2, color=colors[o], label = false)
            
            # Deuxième partie, à partir du point de changement
            plot!(x_line[change_idx:end], vect_line[change_idx:end],
                linestyle=presence_para[change_idx] == 0.0 ? :dashdotdot : :solid, linewidth=2, color=colors[o], label = false)
        end
        #plot!(plt, x_line, vect_line, linewidth=2, color=colors[o], label = labels[o])
    end

    if para_plot_panel[1] == "panel"
        plot!(legend=false)
        plot!(plot_titlefontsize=1)
        plot!(xguidefontsize=10)
        plot!(yguidefontsize=10)

        if nb_modele == 0
            plot!(xlabel = "")
        elseif nb_modele == 1
            plot!(xlabel = "")
            plot!(ylabel = "")
        elseif nb_modele == 2

        elseif nb_modele == 3
            plot!(ylabel = "")
        end
    end 


    #Fixer les paras de sauvgarde
    heure = string(now())
    titre_modele = ["sansPara","ParaOeufs","ParaLarve","ParaNymphe"][nb_modele+1]

    savefig(plt, "$(@__DIR__)/graphs/temporaire/temp.png")
    if YN_save == true
        savefig(plt, joinpath(BIFURCATION_GRAPH_DIR, "Bifurcation_Finale_"*titre_modele*".png"))
        savefig(plt, joinpath(BIFURCATION_GRAPH_DIR, "Bifurcation_Finale_"*titre_modele*".svg"))
    end
    display(plt)

end


### Création de la fonction Affichant le graphique temporaire ----

function printTemp(OnOff=true)
    if OnOff == true
        if Sys.KERNEL == :NT
            run(`cmd /c start "$(@__DIR__)/graphs/temporaire/temp.png"`)       
        else
            run(`open "$(@__DIR__)/graphs/temporaire/temp.png"`)
        end
    end 
end

function Graph_survie_para_1D(bifurction_rez_tot, para_bifucation, vertical,legend, YN_save)

    # Extraire les paramètre pour la figure 
    labels = para_bifucation[5]
    axe_range = exp.(LinRange(log(para_bifucation[2][1]), log(para_bifucation[2][2]), para_bifucation[3]))

    # Les couleurs 
    palette = ColorSchemes.colorschemes[:twelvebitrainbow];
    positions = collect(range(0, 1; length=12));
    indices = [4, 6, 10 ,2];
    indices = indices[1:length(para_bifucation[5])]
    colors = [RGB{Float64}(get(palette, positions[i])) for i in indices];

    α = 0.7  # quantité d'éclaircissement (30%)

    # Évlaircissment des couleurs 
    lighter_colors = [
        RGB(
            c.r + (1 - c.r)*α,
            c.g + (1 - c.g)*α,
            c.b + (1 - c.b)*α
        )
        for c in colors
    ]
    total_colors = vcat(colors, lighter_colors)
    println(length(total_colors))

    # Epaisseur des traits 
    fin = 5
    large = 15
    legende_with = 3

    # Création du plot 
    plt = Plots.plot(
        xlabel="γ : Competition intensity",
        ylabel="", 
        yaxis=false,
        xscale = :log10,
        ylims = (1,length(bifurction_rez_tot[1]) * 3 + 0.5),
        legend=:outerbottom, 
        legendcolumns=2,
        contour_labels=false,
        grid=false);

    # Ajouter les courbes sur le plot
    for o in 1:length(bifurction_rez_tot[1])
        for r in [1,2]
            presence_para = [v[4] for v in bifurction_rez_tot[r][o]] 
            vect_line = fill(o*3-(r-1)*1,para_bifucation[3])
            indice_color = o + length(bifurction_rez_tot[1]) * (r-1)

            # Première valeur différente du vecteur
            change_idx = findfirst(presence_para .!= presence_para[1])

            # Pas de changement
            if isnothing(change_idx)
                Plots.plot!(axe_range, vect_line, 
                    linestyle = presence_para[1] == 0.0 ? :dashdotdot : :solid , 
                    linewidth = presence_para[1] == 0.0 ? fin : large , 
                    color=total_colors[indice_color], label = false)
            # Changement
            else
                # On corrige change_idx pour qu'apres la dernière donné de presence il y ai absance et non apres la première donnée d'absence
                # Première partie, jusqu'au point juste avant le changement
                Plots.plot!(axe_range[1:change_idx-1], vect_line[1:change_idx-1],
                    linestyle = presence_para[1] == 0.0 ? :dashdotdot : :solid, 
                    linewidth = presence_para[1] == 0.0 ? fin : large, 
                    color=total_colors[indice_color], label = false)
                        
                # Deuxième partie, à partir du point de changement
                Plots.plot!(axe_range[change_idx-1:end], vect_line[change_idx-1:end],
                    linestyle = presence_para[change_idx] == 0.0 ? :dashdotdot : :solid, 
                    linewidth = presence_para[change_idx] == 0.0 ? fin : large, 
                    color=total_colors[indice_color], label = false)
                
                if vertical == true 
                    #vline!([mean([axe_range[change_idx-1],axe_range[change_idx]])], color=colors[o], linestyle=:solid, linewidth=1)
                    vline!([axe_range[change_idx]], color=total_colors[indice_color], linestyle=:solid, linewidth=1)
                end
            end
        end
    end 

    # Creation des légendes 
    
    if legend=="zelda" || legend==true  
        for o in 1:length(bifurction_rez_tot[1])       
            # Legend of the colors of the curve
            Plots.plot!(0.00001:0.00001, [v[1] for v in bifurction_rez_tot[1][1]],
                linestyle = :solid, 
                linewidth = legende_with, 
                color=total_colors[o], label = labels[o])
        end 
    end 

    Plots.plot!(0.00001:0.00001, [v[1] for v in bifurction_rez_tot[1][1]], 
        linestyle = :solid , 
        linewidth = legende_with, 
        color=[:black], label = "Hight pest fecundity");

    Plots.plot!(0.00001:0.00001, [v[1] for v in bifurction_rez_tot[1][1]], 
        linestyle = :solid , 
        linewidth = legende_with, 
        color=[:grey], label = "Low pest fecundity");

    Plots.plot!(0.00001:0.00001, [v[1] for v in bifurction_rez_tot[1][1]] , 
            linestyle = :solid , 
            linewidth = legende_with, 
            color=[:black], label = "Survival of the parasitoid");

    Plots.plot!(0.00001:0.00001, [v[1] for v in bifurction_rez_tot[1][1]], 
            linestyle = :dashdotdot , 
            linewidth = legende_with, 
            color=[:black], label = "Death of the parasitoid");

    savefig(plt, "$(@__DIR__)/graphs/temporaire/temp.png")
    if YN_save == true
        savefig(plt, joinpath(FIGURE1D_GRAPH_DIR, "Bifurcation_Finale.png"))
        savefig(plt, joinpath(FIGURE1D_GRAPH_DIR, "Bifurcation_Finale.svg"))
    end
            

    display(plt) 

end


### Création du modèle de dynamique ravageur/parasitoïdes des oeufs ----

# Paramètres 
#   du  : vecteur de float : dérivé de la densité des différentes populations par rapport au temps
#        
#   u   : vecteur de float : densité des différentes populations
#
#   h   : fonction : fonction historique des densité de population
#    
#   p   : vecteur de float : q = fécondité des ravageurs, a = taux d'attaque des parasitoids, alphaL = intensité de la compétition larvaire, 
#           delta = taux de mortalité (U = des larves, L = des oeufs, M = des nymphes, A = des adultes, P = des parasitoïdes),
#           tau = durée du stade (U = des larves, L = des oeufs, M = des nymphes, J = des parasitoïdes juveniles)
#           H et mi non implémenté ici
#
#   t   : float : temps actuel

function without_parasitoid_model(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ, H, mi = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - (tauU + tauL + tauM))[4] 
    hist4 = h(p, t - tauM)[5]
    hist5 = h(p, t - tauL)[2]

    du[1] = q * u[4] - deltaU * u[1] - q * hist1 * exp(- deltaU * tauU)
    du[2] = q * hist1 * exp(- deltaU * tauU)   -   
        q * hist2 * exp(-deltaU * tauU) * exp(u[5])   -   (deltaL + alphaL * u[2]) * u[2] 
    du[3] = q * hist2 * exp(-deltaU * tauU) * exp(u[5]) - deltaM * u[3] - 
        q * hist3 * exp(- deltaU * tauU) * exp(-deltaM * tauM) * exp(hist4)
    du[4] = q * hist3 * exp(- deltaU * tauU) * exp(-deltaM * tauM) * exp(hist4)   -   deltaA * u[4] 
    du[5] = - alphaL * (u[2] - hist5) 
end




### Création du modèle de dynamique ravageur/parasitoïdes des oeufs ----

# Paramètres 
#   du  : vecteur de float : dérivé de la densité des différentes populations par rapport au temps
#        
#   u   : vecteur de float : densité des différentes populations
#
#   h   : fonction : fonction historique des densité de population
#    
#   p   : vecteur de float : q = fécondité des ravageurs, a = taux d'attaque des parasitoids, alphaL = intensité de la compétition larvaire, 
#           delta = taux de mortalité (U = des larves, L = des oeufs, M = des nymphes, A = des adultes, P = des parasitoïdes),
#           tau = durée du stade (U = des larves, L = des oeufs, M = des nymphes, J = des parasitoïdes juveniles)
#           H et mi non implémenté ici
#
#   t   : float : temps actuel

function eggs_parasitoid_model(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ, H, mi = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - tauL)[7] 
    hist4 = h(p, t - (tauU + tauL + tauM))[4]
    hist5 = h(p, t - (tauL + tauM))[7]
    hist6 = h(p, t - tauM)[8]
    hist7 = h(p, t - tauJ)[1]
    hist8 = h(p, t - tauJ)[6]
    hist9 = h(p, t - tauU)[6]
    hist10 = h(p, t - tauL)[2]


    du[1] = q * u[4] - a * u[1] * u[6] - deltaU * u[1] - q * hist1 * exp(u[7])
    du[2] = q * hist1 * exp(u[7]) - (alphaL * u[2] + deltaL) * u[2] - q * hist2 * exp(hist3) * exp(u[8])
    du[3] = q * hist2 * exp(hist3) * exp(u[8]) - deltaM * u[3] - q * hist4 * exp(hist5) * exp(hist6) * exp(- tauM * deltaM)
    du[4] = q * hist4 * exp(hist5) * exp(hist6) * exp(- tauM * deltaM)   -   deltaA * u[4] 
    du[5] = a * u[1] * u[6] - deltaU * u[5] - a * hist7 * hist8 * exp(- deltaU * tauJ)
    du[6] = a * hist7 * hist8 * exp(- deltaU * tauJ) - deltaP * u[6]
    du[7] = - a * (u[6] - hist9) 
    du[8] = - alphaL * (u[2] - hist10)
end





### Création du modèle de dynamique ravageur/parasitoïdes des larves ----

# Paramètres 
#   du  : vecteur de float : dérivé de la densité des différentes populations par rapport au temps
#        
#   u   : vecteur de float : densité des différentes populations
#
#   h   : fonction : fonction historique des densité de population
#    
#   p   : vecteur de float : q = fécondité des ravageurs, a = taux d'attaque des parasitoids, alphaL = intensité de la compétition larvaire, 
#           delta = taux de mortalité (U = des larves, L = des oeufs, M = des nymphes, A = des adultes, P = des parasitoïdes),
#           tau = durée du stade (U = des larves, L = des oeufs, M = des nymphes, J = des parasitoïdes juveniles)
#           H et mi non implémenté ici
#
#   t   : float : temps actuel

function larval_parasitoid_model(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ, H, mi = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - (tauU + tauL + tauM))[4]
    hist4 = h(p, t - tauM)[7]
    hist5 = h(p, t - tauJ)[2]
    hist6 = h(p, t - tauJ)[6]
    hist7 = h(p, t - tauL)[6]
    hist8 = h(p, t - tauL)[2]
    hist9 = h(p, t - tauL)[5]
    hist10 = h(p, t - tauJ)[2]
    hist11 = h(p, t - tauJ)[5]


    du[1] = q * u[4] - deltaU * u[1] - q * hist1 * exp(- deltaU * tauU)
    du[2] = q * hist1 * exp(- deltaU * tauU) - a * u[2] * u[6] - 
       (alphaL * (u[2] + u[5]) + deltaL) * u[2] - q * hist2 * exp(- deltaU * tauU) * exp(u[7])
    du[3] = q * hist2 * exp(- deltaU * tauU) * exp(u[7]) - deltaM * u[3] - 
    q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(- tauM * deltaM)
    du[4] = q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(- tauM * deltaM)   -   deltaA * u[4] 
    du[5] = a * u[2] * u[6] - (alphaL * (u[2] + u[5]) + deltaL) * u[5] - a * hist5 * hist6 * exp(u[8])
    du[6] = a * hist5 * hist6 * exp(u[8]) - deltaP * u[6]
    du[7] = - a * (u[6] - hist7) - alphaL * ((u[2] + u[5]) - (hist8 + hist9))
    du[8] = - alphaL * ((u[2] + u[5])-(hist10 + hist11))
end



### Création du modèle de dynamique ravageur/parasitoïdes des nymphes ----

# Paramètres 
#   du  : vecteur de float : dérivé de la densité des différentes populations par rapport au temps
#        
#   u   : vecteur de float : densité des différentes populations
#
#   h   : fonction : fonction historique des densité de population
#    
#   p   : vecteur de float : q = fécondité des ravageurs, a = taux d'attaque des parasitoids, alphaL = intensité de la compétition larvaire, 
#           delta = taux de mortalité (U = des larves, L = des oeufs, M = des nymphes, A = des adultes, P = des parasitoïdes),
#           tau = durée du stade (U = des larves, L = des oeufs, M = des nymphes, J = des parasitoïdes juveniles)
#           H et mi non implémenté ici
#
#   t   : float : temps actuel

function nymphal_parasitoid_model(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ, H, mi = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - (tauU + tauL + tauM))[4]
    hist4 = h(p, t - tauM)[8]
    hist5 = h(p, t - tauJ)[3]
    hist6 = h(p, t - tauJ)[6]
    hist7 = h(p, t - tauM)[6]
    hist8 = h(p, t - tauL)[2]


    du[1] = q * u[4] - deltaU * u[1] - q * hist1 * exp(- deltaU * tauU) 
    du[2] = q * hist1 * exp(- deltaU * tauU) - (alphaL * u[2] + deltaL) * u[2] -
        q * hist2 * exp(- deltaU * tauU) * exp(u[8])
    du[3] = q * hist2 * exp(- deltaU * tauU) * exp(u[8]) - a * u[3] * u[6] - deltaM * u[3] - 
        q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(u[7])
    du[4] = q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(u[7])   -   deltaA * u[4] 
    du[5] = a * u[3] * u[6] - deltaM * u[5] - a * hist5 * hist6 * exp(- deltaM * tauJ)
    du[6] = a * hist5 * hist6 * exp(- deltaM * tauJ) - deltaP * u[6]
    du[7] = - a * (u[6] - hist7) 
    du[8] = - alphaL * (u[2] - hist8)
end

end



 


