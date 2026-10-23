######################################################################
#### Bifurcation du modèe 1 en comparant avec et sans compétition ####
######################################################################

# Charger les packages et les modules 
include("/Users/Quentin/Desktop/Article1/fonction/fc_dynamique.jl")
using .module_analyse
using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using DataFrames, CSV
using Dates


###########################
## Créaction des modèles ##
###########################


function ssDD_model(du, u, h, p, t)
    alphaL, a, q, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - (tauU + tauL + tauM))[4]
    hist4 = h(p, t - tauM)[7]
    hist5 = h(p, t - tauJ)[2]
    hist6 = h(p, t - tauJ)[6]
    hist7 = h(p, t - tauL)[6]


    du[1] = q * u[4] - deltaU * u[1] - q * hist1 * exp(- tauU * deltaU)
    du[2] = q * hist1 * exp(- tauU * deltaU) - a * u[2] * u[6] - deltaL * u[2] - 
        q * hist2 * exp(- tauU * deltaU) * exp(u[7])
    du[3] = q * hist2 * exp(- tauU * deltaU) * exp(u[7]) - deltaM * u[3] - 
        q * hist3 * exp(- tauU * deltaU) * exp(hist4) * exp(- tauM * deltaM)
    du[4] = q * hist3 * exp(- tauU * deltaU) * exp(hist4) * exp(- tauM * deltaM)  -   deltaA * u[4] 
    du[5] = a * u[2] * u[6] - deltaL * u[5] - a * hist5 * hist6 * exp(- deltaL * tauJ)
    du[6] = a * hist5 * hist6 * exp(- deltaL * tauJ) - deltaP * u[6]
    du[7] = - a * (u[6] - hist7) 
end


function DD_model(du, u, h, p, t)
    alphaL, a, q, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

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


###########################################################
## Création du diagramme de bifurcation sans compétition ## 
###########################################################


# Création paramètres et fixation de leur paramètres

# Création des paramètres du modèle
DicoDesParasGlobal = Dict(
    :"tauU" => 4.72,
    :"tauL" => 16,
    :"tauM" => 8.92,
    :"tauJ" => 4.72,
    :"alphaL" =>  0.001,
    :"deltaU" =>  0.0432,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"deltaA" => 1/12.2,
    :"deltaP" =>  1,
    :"a" => 1,
    :"q" => 106.4/12.2,
)

# Parametre de simulation
tspan = (0.0, 20000.0)

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("alphaL","a", "q", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para) 
    return [1, 1, 1, 1, 1, 1, - dico_para["tauL"] * (dico_para["a"] + dico_para["deltaL"])]
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
   return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]), 
   (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"]),
   dico_para["tauM"],dico_para["tauJ"],dico_para["tauJ"],dico_para["tauL"]]
end


# para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
para_du = (filter(kv -> kv[1] != "deltaU", DicoDesParasGlobal),fct_lags, ordre_para)
para_dl = (filter(kv -> kv[1] != "deltaL", DicoDesParasGlobal),fct_lags, ordre_para)
para_dm = (filter(kv -> kv[1] != "deltaM", DicoDesParasGlobal),fct_lags, ordre_para)
para_da = (filter(kv -> kv[1] != "deltaA", DicoDesParasGlobal),fct_lags, ordre_para)
    

# para_run 
    # fonction créant les valeurs d'initialisation des variables 
    # début et fin de la dynamique 
    # méthode de résolution numérique
para_run = (fct_hini,tspan,MethodOfSteps(Tsit5()))

# para_matrice
    # tuple des deux noms de paramètres qui serviront d'axe à la matrice
    # tuple des valeur min et max des deux paramètres               
    # tuple des dimensions des deux axes de la matrice
    # int des X dernier pas de temps étudié
para_bifucation_du = ("deltaU",(0,1),30,1000)
para_bifucation_dl = ("deltaL",(0,1),30,1000)
para_bifucation_dm = ("deltaM",(0,1),30,1000)
para_bifucation_da = ("deltaA",(0,1),30,1000)
    
    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (2,false,true) 
    
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # affiche toutes les valeurs prise ou juste les max/min 
para_plot = (1,true)
    
# Appelle de la fonction dynamique 
@time plotUssDD = module_analyse.simulate_bifurcation(
    ssDD_model, para_du, para_run, para_bifucation_du, para_save, para_plot)
@time plotLssDD = module_analyse.simulate_bifurcation(
    ssDD_model, para_dl, para_run, para_bifucation_dl, para_save, para_plot)
@time plotMssDD = module_analyse.simulate_bifurcation(
    ssDD_model, para_dm, para_run, para_bifucation_dm, para_save, para_plot)
#@time plotAssDD = module_analyse.simulate_bifurcation(
#    ssDD_model, para_da, para_run, para_bifucation_da, para_save, para_plot)




###########################################################
## Création du diagramme de bifurcation avec compétition ## 
###########################################################

# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para) 
    return [1, 1, 1, 1, 1, 1, 
       - dico_para["tauL"] * (2 * dico_para["alphaL"] + dico_para["deltaL"] + dico_para["a"]),
       - dico_para["tauJ"] * (2 * dico_para["alphaL"] + dico_para["deltaL"])]
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
   return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]),
   (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"]),
   dico_para["tauM"],dico_para["tauJ"], dico_para["tauJ"], dico_para["tauL"],dico_para["tauL"],
   dico_para["tauJ"],dico_para["tauJ"]]
end

# para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
    para_du = (filter(kv -> kv[1] != "deltaU", DicoDesParasGlobal),fct_lags, ordre_para)
    para_dl = (filter(kv -> kv[1] != "deltaL", DicoDesParasGlobal),fct_lags, ordre_para)
    para_dm = (filter(kv -> kv[1] != "deltaM", DicoDesParasGlobal),fct_lags, ordre_para)
    para_da = (filter(kv -> kv[1] != "deltaA", DicoDesParasGlobal),fct_lags, ordre_para)
        
    
    # para_run 
        # fonction créant les valeurs d'initialisation des variables 
        # début et fin de la dynamique 
        # méthode de résolution numérique
    para_run = (fct_hini,tspan,MethodOfSteps(Tsit5()))

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
    para_save = (2,true,true) 


# Appelle de la fonction dynamique 
@time plotUDD = module_analyse.simulate_bifurcation(
    DD_model, para_du, para_run, para_bifucation_du, para_save, para_plot)
@time plotLDD = module_analyse.simulate_bifurcation(
    DD_model, para_dl, para_run, para_bifucation_dl, para_save, para_plot)
@time plotMDD = module_analyse.simulate_bifurcation(
    DD_model, para_dm, para_run, para_bifucation_dm, para_save, para_plot)
#@time plotADD = module_analyse.simulate_bifurcation(
#    DD_model, para_da, para_run, para_bifucation_da, para_save, para_plot)



########################
## Affichage fusionné ##
########################

#Information générale
axe_range = range(0,1,30)
 yname = "Valeurs asymptotiques"

# Birfurcation pour deltaU
pltU = plot(title="Diagramme de bifurcation Modèle 1 en fonction de deltaU", 
xlabel= "DeltaU", 
ylabel = yname)

for (i, axe) in enumerate(axe_range)
    scatter!(pltU, 
        fill(axe, length(plotUssDD[i])),  # Toutes les valeurs de q (fixe) pour les résultats correspondants
        plotUssDD[i],                   # Les résultats (valeurs asymptotiques ou amplitudes)
        label=true,                  
        color=:black,                 
        alpha=0.003,                  
        ms=2)      
     scatter!(pltU, 
        fill(axe, length(plotUDD[i])),  # Toutes les valeurs de q (fixe) pour les résultats correspondants
        plotUDD[i],                   # Les résultats (valeurs asymptotiques ou amplitudes)
        label=true,                 
        color=:red,                 
        alpha=0.003,                  
        ms=2)                    
end

display(pltU)

function doubleBifurcation(dimension, nbmodele, para_bif, data1, data2)
    axe_range = range(dimension)
    yname = "Valeurs asymptotiques"
    
    pltU = plot(title="Diagramme de bifurcation Modèle "*string(nbmodele)*" en fonction de "*string(para_bif)*" ", 
    xlabel= string(para_bif), 
    ylabel = yname)

    for (i, axe) in enumerate(axe_range)
        scatter!(pltU, 
            fill(axe, length(data1[i])),  # Toutes les valeurs de q (fixe) pour les résultats correspondants
            data1[i],                   # Les résultats (valeurs asymptotiques ou amplitudes)
            label=true,                  
            color=:black,                 
            alpha=0.003,                  
            ms=2)      
        scatter!(pltU, 
            fill(axe, length(data2[i])),  # Toutes les valeurs de q (fixe) pour les résultats correspondants
            data2[i],                   # Les résultats (valeurs asymptotiques ou amplitudes)
            label=true,                 
            color=:red,                 
            alpha=0.003,                  
            ms=2)                    
    end
end

# Fusionner les bifurcations 
heure = string(now())

plotU = plot(plotUssDD)
plot!(plotU, plotUDD)
savefig(plt, "/Users/Quentin/Desktop/Article1/graphs/ParaOeufs/dataBifurcation_deltaU_Total_"*heure*".png")      

plotL = plot(plotLssDD)
plot!(plotL, plotLDD)
savefig(plt, "/Users/Quentin/Desktop/Article1/graphs/ParaOeufs/dataBifurcation_deltaL_Total_"*heure*".png")      


plotM = plot(plotMssDD)
plot!(plotM, plotMDD)
savefig(plt, "/Users/Quentin/Desktop/Article1/graphs/ParaOeufs/dataBifurcation_deltaM_Total_"*heure*".png")      

