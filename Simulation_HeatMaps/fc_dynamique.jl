######################################################################################################
#### Répertoire des fonction pour la visualisation et l'analyse des fonction des modèles à retard ####
######################################################################################################

### Code générale ---- 
module module_analyse

# Fonction exporté 
export  heatMapStab, printTemp, dicoPara,IntegrerARGS, fixeAxes


# Import packages
using CSV, Trapz, FFTW,DifferentialEquations
using DataFrames, Statistics, Dates


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
    
    # Commencer le timer
    start_time = time()


    # Extraire les paramètres du modèle 
    fixed_params = para[1]
    lags_fct =  para[2]
    ordre_para = para[3]

    # Extraire les paramètres de simulation
    tspan = para_run[2]
    alg = para_run[3]
    pdt = para_run[4]

    t_stud = para_matrice[4]

    # Extraire les paramètres de sauvgarde
    nb_modele = para_save[1]
    presence_dd = para_save[2]
    YN_save = para_save[3]
    SaveHisto = para_save[4]
    petit_nom = para_save[5]
    if petit_nom == "pdt"
        petit_nom = string(pdt)
    end

    # Extraire les paramètres pour le graphique
    pop_obs = para_plot[1]
    YN_stab = para_plot[2]
    YN_mean = para_plot[3]
    

    # Etraire les paramètre pour la création de la heapmap
    axe1_name = string(para_matrice[1][1])
    axe2_name = string(para_matrice[1][2])
    axe1_range = range(para_matrice[2][1][1],para_matrice[2][1][2],para_matrice[3][1])
    axe2_range = range(para_matrice[2][2][1],para_matrice[2][2][2],para_matrice[3][2])

    resultsVar = zeros(para_matrice[3][1], para_matrice[3][2])
    resultMean = zeros(para_matrice[3][1], para_matrice[3][2])
    resultOeuf = zeros(para_matrice[3][1], para_matrice[3][2])
    resultNymphe = zeros(para_matrice[3][1], para_matrice[3][2])
    resultAdulte = zeros(para_matrice[3][1], para_matrice[3][2])
    resultParaJuv = zeros(para_matrice[3][1], para_matrice[3][2])
    resultPara = zeros(para_matrice[3][1], para_matrice[3][2])
    resultsMax = zeros(para_matrice[3][1], para_matrice[3][2])
    resultsMin = zeros(para_matrice[3][1], para_matrice[3][2])
    resultsFreq = zeros(para_matrice[3][1], para_matrice[3][2])
    resultsIntegral = zeros(para_matrice[3][1], para_matrice[3][2])
    SaveHisto && (resultsHisto = fill(Vector{Float64}(), para_matrice[3][1] * para_matrice[3][2]))




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
                        recent_values = [v[2] for v in integrator.sol.u[end-(length(recent_time)-1):end]]
                        std_l = std(recent_values)

                        if std_l < 0.0001
                            return true
                        elseif (abs(maximum(recent_values[Int(round(length(recent_values)/2)):end]) - 
                            maximum(recent_values[1:Int(round(length(recent_values)/2))-1]))) / mean(recent_values) < 0.0001 &&
                            abs(minimum(recent_values[Int(round(length(recent_values)/2)):end]) - 
                            minimum(recent_values[1:Int(round(length(recent_values)/2))-1])) / mean(recent_values) < 0.0001
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
            elseif alg == "fixeSimple"
                # Définir le problème et le résoudre
                frq_save = 0.1
                prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
                sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=frq_save,abstol=1e-9, reltol=1e-9) 

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

            # Sélectionner les dernières unités de temps 
            t_values = collect(t_values)

            if t_values[end] > t_stud
                t_last_values = [x for x in t_values if x > (t_values[end] - t_stud)]
                u_last_values = u_values[end-length(t_last_values)+1:end] 
            else 
                t_last_values = t_values
                u_last_values = u_values
            end


            # Selectionner la population observé
            if length(pop_obs)<2
                u_obs = [v[pop_obs] for v in u_last_values]
                resultMean[j, i] = mean(u_obs)
                resultOeuf[j, i] = 0
                resultNymphe[j, i] = 0
                resultAdulte[j, i] = 0
                resultParaJuv[j, i] = 0
                resultPara[j, i] = 0

            elseif length(pop_obs)==6
                o_obs = [v[pop_obs[1]] for v in u_last_values]
                u_obs = [v[pop_obs[2]] for v in u_last_values]
                n_obs = [v[pop_obs[3]] for v in u_last_values]
                a_obs = [v[pop_obs[4]] for v in u_last_values]
                j_obs = [v[pop_obs[5]] for v in u_last_values]
                p_obs = [v[pop_obs[6]] for v in u_last_values]
                resultMean[j, i] = mean(u_obs)
                resultOeuf[j, i] = mean(o_obs)
                resultNymphe[j, i] = mean(n_obs)
                resultAdulte[j, i] = mean(a_obs)
                resultParaJuv[j, i] = mean(j_obs)
                resultPara[j, i] = mean(p_obs)

            elseif length(pop_obs)==4
                o_obs = [v[pop_obs[1]] for v in u_last_values]
                u_obs = [v[pop_obs[2]] for v in u_last_values]
                n_obs = [v[pop_obs[3]] for v in u_last_values]
                a_obs = [v[pop_obs[4]] for v in u_last_values]
                resultMean[j, i] = mean(u_obs)
                resultOeuf[j, i] = mean(o_obs)
                resultNymphe[j, i] = mean(n_obs)
                resultAdulte[j, i] = mean(a_obs)
                resultParaJuv[j, i] = 0
                resultPara[j, i] = 0
            end

            # Calculer la variance et la moyenne des 100 derniers pas de temps

            resultsVar[j, i] = std(u_obs)
            resultsMax[j, i] = maximum(u_obs)
            resultsMin[j, i] = minimum(u_obs)
            SaveHisto && (resultsHisto[Int(round(para_matrice[3][1]*(i-1)+j))] = 
                [u_obs[i] for i in 1:length(t_last_values) if t_last_values[i]%1==0])
            SaveHisto && (append!(resultsHisto[Int(round(para_matrice[3][1]*(i-1)+j))],[Float64(j),Float64(i)])) 
            
            if resultsVar[j, i] > 0.0001
                # Trouver les indices des pics locaux
                peak_indices = findall(i -> u_obs[i-1] < u_obs[i] > u_obs[i+1], 2:length(u_obs)-1) .+ 1  # Corriger les indices

                # Calculer les écarts de temps entre les pics
                if length(peak_indices) > 1
                    periods = diff(t_last_values[peak_indices])  # Différence entre les instants des pics
                    dominant_period = mean(periods)  # Moyenne des périodes détectées
                    resultsFreq[j, i] = dominant_period  
                    
                    # Trouver resultsIntegral
                    t0 = t_last_values[end] - dominant_period  # Prendre la dernière période pour être sûr d’avoir un cycle complet
                    indices = findall(x -> t0 ≤ x ≤ t_last_values[end], t_last_values)
                    resultsIntegral[j, i] = trapz(t_last_values[indices], u_obs[indices])
                end 
            end
        end
    end

    # Nettoyage des données invalides
    resultsVar[resultMean .== 0] .= 0

    if YN_save == true
         
        # Fixer les paras de sauvgarde
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
        resultOeuf = DataFrame(resultPara, Symbol.(axe1_range))
        resultNymphe = DataFrame(resultPara, Symbol.(axe1_range))
        resultAdulte = DataFrame(resultPara, Symbol.(axe1_range))
        resultParaJuv = DataFrame(resultPara, Symbol.(axe1_range))
        resultPara = DataFrame(resultPara, Symbol.(axe1_range))
        resultsMax = DataFrame(resultsMax, Symbol.(axe1_range))
        resultsMin = DataFrame(resultsMin, Symbol.(axe1_range))        
        resultsFreq = DataFrame(resultsFreq, Symbol.(axe1_range))
        resultsIntegral = DataFrame(resultsIntegral, Symbol.(axe1_range))

        # Ajouter une colonne pour les légendes des lignes
        resultsVar[!, "Row Labels"] = axe2_range
        resultMean[!, "Row Labels"] = axe2_range
        resultOeuf[!, "Row Labels"] = axe2_range
        resultNymphe[!, "Row Labels"] = axe2_range
        resultAdulte[!, "Row Labels"] = axe2_range
        resultParaJuv[!, "Row Labels"] = axe2_range
        resultPara[!, "Row Labels"] = axe2_range

        nom_modele = ["sansPara","ParaOeufs","ParaLarve","ParaNymphe"][nb_modele+1]
        cd(@__DIR__)

        # Créer un dossier
        dossier_path = "REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"
        if !isdir(dossier_path)
            mkpath(dossier_path)
        end

        # Sauvgarder dans le dossier
        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Variance_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsVar) 
        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Densité_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultMean)
        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Densité_Oeuf_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultOeuf)
        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Densité_Nymphe_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultNymphe)
        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Densité_Adulte_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultAdulte)
        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Densité_ParaJuv_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultParaJuv)
        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Densité_Para_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultPara)

        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Frequence_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsFreq;)
        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Integral_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsIntegral;) 
        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Maximum_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsMax;)
        CSV.write("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_heatmap_Minimum_"*DD*"_"*axe1_name*"_"*axe2_name*".csv", resultsMin;)
        
        cd(@__DIR__)
        if SaveHisto 
            open("REZ/"*petit_nom*"/"*string(alg)*"/"*string(pdt)*"/"*nom_modele*"_valeur_historique_"*DD*"_"*axe1_name*"_"*axe2_name*".txt", "w") do file
                for subvector in resultsHisto
                    write(file, join(subvector, ",") * "\n")
                end
            end
        end

    end
    
    # Finir le timer et sortir le temps d'exécution
    end_time = time()
    println("Durée d'exécution : ", end_time - start_time, " s")


    return resultMean, resultsVar 
end

function printTemp()
    cd(@__DIR__)
    if Sys.KERNEL == :NT
        run(`cmd /c start "REZ/temporaire/temp.png"`)       
    else
        run(`open "REZ/temporaire/temp.png"`)
    end 
end

function dicoPara(type="base")
    if type == "base"
        Dico = Dict(
            :"deltaU" =>  0.04,
            :"tauU" => 5,
            :"tauL" => 16,
            :"tauM" => 9,
            :"deltaL" => 0.005,
            :"deltaM" => 0.02,
            :"deltaA" => 0.1,
            :"deltaP" =>  0.1,
            :"a" => 0.01,
            :"alphaL" => 0.001,
            :"q" => 9,)
    elseif type == "nul"
        Dico = Dict(
            :"deltaU" =>  0.0,
            :"tauU" => 4,
            :"tauL" => 16,
            :"tauM" => 8,
            :"deltaL" => 0.0,
            :"deltaM" => 0.0,
            :"deltaA" => 0.1,
            :"deltaP" =>  0.1,
            :"a" => 0.01,
            :"alphaL" => 0.0001,
            :"q" => 8,)
    elseif type == "egale"
        Dico = Dict(
            :"deltaU" =>  0.0,
            :"tauU" => 4,
            :"tauL" => 4,
            :"tauM" => 4,
            :"deltaL" => 0.0,
            :"deltaM" => 0.0,
            :"deltaA" => 0.1,
            :"deltaP" =>  0.1,
            :"a" => 0.01,
            :"alphaL" => 0.0001,
            :"q" => 8,)
    end 

    return Dico
end

function IntegrerARGS(arg)

    # Vérification du fonctionement du array 
    # Vérifier que les arguments sont fournis
    if length(arg) < 6
        error("Veuillez fournir deux arguments.")
    end

    # Convertir les arguments en nombres
    numSimul = parse(Int64, arg[1])
    nbModeles = parse(Int64, arg[2])
    pdt = parse(Float64, arg[3])
    alg = arg[4]
    dimHM = parse(Int64, arg[5])
    if parse(Int64, arg[6])==1 
        histSave = true
    elseif parse(Int64, arg[6])==0
        histSave = false
    end
    ptnom = arg[7]
    typeSimul = (numSimul ÷ nbModeles) + 1 
    if (numSimul % nbModeles == 0) typeSimul = (numSimul ÷ nbModeles) end 

    # Paramètres de simulation
    tspan = (0.0, 20000.0)

    return typeSimul, pdt, alg, dimHM, tspan, histSave, ptnom
end

function fixeAxes(nbtypes,fixAbs, fixOrd,Dico)
    # Initialiser un vecteur vide pour stocker les dictionnaires
    dict_vector = Vector{Dict{String, Any}}()

    # Boucle création des dico
    # TODO

end

end