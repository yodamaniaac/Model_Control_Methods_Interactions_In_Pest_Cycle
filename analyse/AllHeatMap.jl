##################################################
#### Analyse Modèle Ravageur tout les modèles ####
##################################################


# Charger les packages et les modules 

if Sys.KERNEL == :NT
    include("P:/Article1/fonction/fc_Solveur.jl")
    include("P:/Article1/fonction/fc_dynamique.jl")
    include("P:/Article1/fonction/fc_plots.jl")
else
    include("/Users/Quentin/Desktop/Article1/fonction/fc_dynamique.jl")
    include("/Users/Quentin/Desktop/Article1/fonction/fc_Solveur.jl")
    include("/Users/Quentin/Desktop/Article1/fonction/fc_plots.jl")
end

using .solveur_maison
using .module_analyse
using .module_graphique
using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using DataFrames, CSV
using Dates
using Statistics





#########################################
## Créaction de l'ensemble des modèles ##
#########################################

function model_ssPara(du, u, h, p, t)
    q, alphaL, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM, tau1, tau2, tau3, tau4, tau5 = p

    hist1 = h(p, t - tau1)[2]
    hist2 = h(p, t - tau2)[2]
    hist3 = h(p, t - tau3)[2] 
    hist4 = h(p, t - tau4)[3]
    hist5 = h(p, t - tau5)[1]

    du[1] = q * hist1 * exp(- deltaU * tauU)   -   
        q * hist2 * exp(-deltaU * tauU) * exp(u[3])   -   (deltaL + alphaL * u[1]) * u[1] 
    du[2] = q * hist3 * exp(- deltaU * tauU) * exp(-deltaM * tauM) * exp(hist4)   -   deltaA * u[2] 
    du[3] = - alphaL * (u[1] - hist5) 
end

function model1_DD(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

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
    du[3] = q * hist2 * exp(hist3) * exp(u[8]) - deltaM * u[3] - 
    q * hist4 * exp(hist5) * exp(hist6) * exp(- tauM * deltaM)
    du[4] = q * hist4 * exp(hist5) * exp(hist6) * exp(- tauM * deltaM)   -   deltaA * u[4] 
    du[5] = a * u[1] * u[6] - deltaU * u[5] - a * hist7 * hist8 * exp(- deltaU * tauJ)
    du[6] = a * hist7 * hist8 * exp(- deltaU * tauJ) - deltaP * u[6]
    du[7] = - a * (u[6] - hist9) 
    du[8] = - alphaL * (u[2] - hist10)
end

function model1_ssDD(du, u, h, p, t)
    q, a, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - tauL)[7] 
    hist4 = h(p, t - (tauU + tauL + tauM))[4]
    hist5 = h(p, t - (tauL + tauM))[7]
    hist6 = h(p, t - tauJ)[1]
    hist7 = h(p, t - tauJ)[6]
    hist8 = h(p, t - tauU)[6]


    du[1] = q * u[4] - a * u[1] * u[6] - deltaU * u[1] - q * hist1 * exp(u[7])
    du[2] = q * hist1 * exp(u[7]) - deltaL * u[2] - q * hist2 * exp(hist3) * exp(- tauL * deltaL)
    du[3] = q * hist2 * exp(hist3) * exp(- tauL * deltaL) - deltaM * u[3] - q * hist4 * exp(hist5) * exp(- tauL * deltaL) * exp(- tauM * deltaM)
    du[4] = q * hist4 * exp(hist5) * exp(- tauL * deltaL) * exp(- tauM * deltaM)   -   deltaA * u[4] 
    du[5] = a * u[1] * u[6] - deltaU * u[5] - a * hist6 * hist7 * exp(- deltaU * tauJ)
    du[6] = a * hist6 * hist7 * exp(- deltaU * tauJ) - deltaP * u[6]
    du[7] = - a * (u[6] - hist8) 
end




function model2_DD(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p
    

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


function model2_ssDD(du, u, h, p, t)
    q, a, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p
    

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



function model3_DD(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p
    

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





function model3_ssDD(du, u, h, p, t)
    q, a, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p
    

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - (tauU + tauL + tauM))[4]
    hist4 = h(p, t - tauJ)[3]
    hist5 = h(p, t - tauJ)[6]
    hist6 = h(p, t - tauM)[6]


    du[1] = q * u[4] - deltaU * u[1] - q * hist1 * exp(- tauU * deltaU)
    du[2] = q * hist1 * exp(- tauU * deltaU) - deltaL * u[2] - 
        q * hist2 * exp(- tauU * deltaU) * exp(- tauL * deltaL)
    du[3] = q * hist2 * exp(- tauU * deltaU) * exp(- tauL * deltaL) - a * u[3] * u[6] - deltaM * u[3] - 
        q * hist3 * exp(- tauU * deltaU)* exp(- tauL * deltaL) * exp(u[7]) 
    du[4] = q * hist3 * exp(- tauU * deltaU)* exp(- tauL * deltaL) * exp(u[7])  -   deltaA * u[4] 
    du[5] = a * u[3] * u[6] - deltaM * u[5] - a * hist4 * hist5 * exp(- deltaM * tauJ)
    du[6] = a * hist4 * hist5 * exp(- deltaM * tauJ) - deltaP * u[6]
    du[7] = - a * (u[6] - hist6) 
end


############################################
## Tracer les heatmap de tout les modèles ##
############################################


# Modèles sans Parasitoïdes

# Création des paramètres du modèle

# Création du Dico global
DicoDesParasGlobal = Dict(
    :"tauU" => 4.72,
    :"tauL" => 16,
    :"tauM" => 8.92,
    :"tau1" =>  4.72,
    :"tau2" => 4.72 + 16,
    :"tau3" =>  4.72 + 16 + 8.92,
    :"tau4" =>  8.92,
    :"tau5" =>  16,
    :"alphaL" =>  0.001,
    :"deltaA" => 1/12.2,
    :"deltaU" =>  0.0432,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"q" => 106.4/12.2,
)


# Paramètres de simulation
tspan = (0.0, 500000.0)
alg = "maison"
pdt = 0.01

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "tauU", "tauL", "tauM", "tau1", "tau2", "tau3", "tau4", "tau5")


# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        return [0, 0, - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
    else 
        return [1, 1, - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
    end 
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
   return [dico_para["tau1"],dico_para["tau2"],dico_para["tau3"],dico_para["tau4"],dico_para["tau5"]]
end


# para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
para_du = (filter(kv -> kv[1] !="deltaU", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)
para_dl = (filter(kv -> kv[1] !="deltaL", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)
para_dm = (filter(kv -> kv[1] !="deltaM", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)
para_da = (filter(kv -> kv[1] !="deltaA", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)

# para_run 
    # fonction créant les valeurs d'initialisation des variables 
    # début et fin de la dynamique 
    # méthode de résolution numérique
para_run = (fct_hini,tspan,alg,pdt)

# para_matrice
    # tuple des deux noms de paramètres qui serviront d'axe à la matrice
    # tuple des valeur min et max des deux paramètres               
    # tuple des dimensions des deux axes de la matrice
    # int des X dernier pas de temps étudié
para_matrice_du = (("q","deltaU"),((1,70),(0.001,1)),(5,5),200)
para_matrice_dl = (("q","deltaL"),((1,70),(0.001,1)),(5,5),200)
para_matrice_dm = (("q","deltaM"),((1,70),(0.001,1)),(5,5),200)
para_matrice_da = (("q","deltaA"),((1,70),(0.001,1)),(5,5),200)

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (0,true,true) 
    
        
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # tracer la heatmap de stabilité
        # tracer la heatmap des moyennes 
para_plot = (1,true, true)


# Appelle de la fonction dynamique 
@time module_analyse.heatMapStab(
    model_ssPara, para_du, para_run, para_matrice_du, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model_ssPara, para_dl, para_run, para_matrice_dl, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model_ssPara, para_dm, para_run, para_matrice_dm, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model_ssPara, para_da, para_run, para_matrice_da, para_save, para_plot)
    module_analyse.printTemp()












# Modèle Parasitoïde des oeufs avec compétition

# Création des paramètres du modèle

# Création du Dico global
DicoDesParasGlobal = Dict(
    :"deltaU" =>  0.0432,
    :"tauU" => 4.72,
    :"tauL" => 16,
    :"tauM" => 8.92,
    :"tauJ" => 4.72,
    :"alphaL" =>  0.001,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"deltaA" => 1/12.2,
    :"deltaP" =>  1,
    :"a" => 1,
    :"alphaL" => 0.001,
    :"q" => 106.4/12.2,
)


# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q","a", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")

# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        return [1, 1, 1, 1, 1, 1, - dico_para["tauU"] * (dico_para["a"] + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
    else 
        return [0, 0, 0, 0, 0, 0, - dico_para["tauU"] * (dico_para["a"] + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
    end 
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
    return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]), dico_para["tauL"],
    (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"]),
    (dico_para["tauL"] + dico_para["tauM"]),dico_para["tauM"],dico_para["tauJ"],
    dico_para["tauJ"],dico_para["tauU"],dico_para["tauL"]]
end

        
    # para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
    para_du = (filter(kv -> kv[1] !="deltaU", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)
    para_dl = (filter(kv -> kv[1] !="deltaL", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)
    para_dm = (filter(kv -> kv[1] !="deltaM", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)
    para_da = (filter(kv -> kv[1] !="deltaA", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
        para_save = (1,true,true) 

    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # tracer la heatmap de stabilité
        # tracer la heatmap des moyennes 
para_plot = (2,true, true)


# Appelle de la fonction dynamique 
@time module_analyse.heatMapStab(
    model1_DD, para_du, para_run, para_matrice_du, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model1_DD, para_dl, para_run, para_matrice_dl, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model1_DD, para_dm, para_run, para_matrice_dm, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model1_DD, para_da, para_run, para_matrice_da, para_save, para_plot)
    module_analyse.printTemp()












# Modèle Parasitoïde des oeufs sans compétition

# Création des paramètres du modèle

# Création du Dico global
DicoDesParasGlobal = Dict(
    :"deltaU" =>  0.0432,
    :"tauU" => 4.72,
    :"tauL" => 16,
    :"tauM" => 8.92,
    :"tauJ" => 4.72,
    :"alphaL" =>  0.001,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"deltaA" => 1/12.2,
    :"deltaP" =>  1,
    :"a" => 1,
    :"alphaL" => 0.001,
)


# Structurer les paramètres pour les inserer dans le modèle
 ordre_para = ("q","a", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")

# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        return [1, 1, 1, 1, 1, 1, - dico_para["tauU"] * (dico_para["a"] + dico_para["deltaU"])]
    else 
        return [0, 0, 0, 0, 0, 0, - dico_para["tauU"] * (dico_para["a"] + dico_para["deltaU"])]
    end 
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
    return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]), dico_para["tauL"],
    (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauL"]),
    (dico_para["tauL"] + dico_para["tauL"]),dico_para["tauJ"],dico_para["tauJ"],dico_para["tauU"]]
 end

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
        para_save = (1,false,true) 

# Appelle de la fonction dynamique 
@time module_analyse.heatMapStab(
    model1_ssDD, para_du, para_run, para_matrice_du, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model1_ssDD, para_dl, para_run, para_matrice_dl, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model1_ssDD, para_dm, para_run, para_matrice_dm, para_save, para_plot)
    module_analyse.printTemp() 


@time module_analyse.heatMapStab(
    model1_ssDD, para_da, para_run, para_matrice_da, para_save, para_plot)
    module_analyse.printTemp()










# Modèle Parasitoïde des larves avec compétition

# Création des paramètres du modèle

# Création du Dico global
DicoDesParasGlobal = Dict(
    :"tauU" => 4.72,
    :"tauL" => 16,
    :"tauM" => 8.92,
    :"tauJ" => 16,
    :"alphaL" =>  0.001,
    :"deltaA" => 1/12.2,
    :"deltaU" =>  0.0432,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"q" => 106.4/12.2,
    :"a" => 1,
)



# Structurer les paramètres pour les inserer dans le modèle
 ordre_para = ("q","a", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")

# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        return [1, 1, 1, 1, 1, 1, 
        - dico_para["tauL"] * (2 * dico_para["alphaL"] + dico_para["deltaL"] + dico_para["a"]),
        - dico_para["tauJ"] * (2 * dico_para["alphaL"] + dico_para["deltaL"])]
    else 
        return [0, 0, 0, 0, 0, 0, 
        - dico_para["tauL"] * (2 * dico_para["alphaL"] + dico_para["deltaL"] + dico_para["a"]),
        - dico_para["tauJ"] * (2 * dico_para["alphaL"] + dico_para["deltaL"])]
    end 
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
    return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]),
    (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"]),
    dico_para["tauM"],dico_para["tauJ"], dico_para["tauJ"], dico_para["tauL"],dico_para["tauL"],
    dico_para["tauJ"],dico_para["tauJ"]]
end

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
        para_save = (2,true,true) 

# Appelle de la fonction dynamique 
@time module_analyse.heatMapStab(
    model2_DD, para_du, para_run, para_matrice_du, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model2_DD, para_dl, para_run, para_matrice_dl, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model2_DD, para_dm, para_run, para_matrice_dm, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model2_DD, para_da, para_run, para_matrice_da, para_save, para_plot)
    module_analyse.printTemp()












# Modèle Parasitoïde des larves avec compétition

# Création des paramètres du modèle

# Création du Dico global
DicoDesParasGlobal = Dict(
    :"q" => 104.6/12/2,
    :"tauU" => 4.72,
    :"tauL" => 16,
    :"tauM" => 8.92,
    :"tauJ" => 16,
    :"alphaL" =>  0.001,
    :"deltaU" => 0.0432,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"deltaA" => 1/12.2,
    :"deltaP" =>  1,
    :"a" => 1,
)

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q","a", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")

# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        return [1, 1, 1, 1, 1, 1, - dico_para["tauL"] * (dico_para["a"] + dico_para["deltaL"])]
    else 
        return [0, 0, 0, 0, 0, 0, - dico_para["tauL"] * (dico_para["a"] + dico_para["deltaL"])]
    end 
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
    return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]), 
    (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"]),
    dico_para["tauM"],dico_para["tauJ"],dico_para["tauJ"],dico_para["tauL"]]
 end


     # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
        para_save = (2,false,true) 

# Appelle de la fonction dynamique 
@time module_analyse.heatMapStab(
    model2_ssDD, para_du, para_run, para_matrice_du, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model2_ssDD, para_dl, para_run, para_matrice_dl, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model2_ssDD, para_dm, para_run, para_matrice_dm, para_save, para_plot)
    module_analyse.printTemp() 


@time module_analyse.heatMapStab(
    model2_ssDD, para_da, para_run, para_matrice_da, para_save, para_plot)
    module_analyse.printTemp()













# Modèle Parasitoïde des nymphes avec compétition

# Création des paramètres du modèle

# Création du Dico global
DicoDesParasGlobal = Dict(
    :"tauU" => 4.72,
    :"tauL" => 16,
    :"tauM" => 8.92,
    :"tauJ" => 8.92,
    :"alphaL" =>  0.001,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"deltaA" => 1/12.2,
    :"deltaP" =>  1,
    :"a" => 1,
    :"q" => 106.4/12.2,
    :"alphaL" => 0.001,
)


# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q","a", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")

# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        return [1, 1, 1, 1, 1, 1, - dico_para["tauU"] * (dico_para["a"] + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
    else 
        return [0, 0, 0, 0, 0, 0, - dico_para["tauU"] * (dico_para["a"] + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
    end 
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
    return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]), dico_para["tauL"],
    (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"]),
    (dico_para["tauL"] + dico_para["tauM"]),dico_para["tauM"],dico_para["tauJ"],
    dico_para["tauJ"],dico_para["tauU"],dico_para["tauL"]]
 end


     # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
        para_save = (3,true,true) 


# Appelle de la fonction dynamique 
@time module_analyse.heatMapStab(
    model3_DD, para_du, para_run, para_matrice_du, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model3_DD, para_dl, para_run, para_matrice_dl, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model3_DD, para_dm, para_run, para_matrice_dm, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model3_DD, para_da, para_run, para_matrice_da, para_save, para_plot)
    module_analyse.printTemp()










# Modèle Parasitoïde des nymphes sans compétition

# Création des paramètres du modèle

# Création du Dico global
DicoDesParasGlobal = Dict(
    :"q" => 104.6/12/2,
    :"tauU" => 4.72,
    :"tauL" => 16,
    :"tauM" => 8.92,
    :"tauJ" => 8.92,
    :"alphaL" =>  0.001,
    :"deltaU" => 0.0432,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"deltaA" => 1/12.2,
    :"deltaP" =>  1,
    :"a" => 1,
)


# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q","a", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")

# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        return [1, 1, 1, 1, 1, 1, - dico_para["tauM"] * (dico_para["a"] + dico_para["deltaM"])]
    else 
        return [0, 0, 0, 0, 0, 0, - dico_para["tauM"] * (dico_para["a"] + dico_para["deltaM"])]
    end 
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
    return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]), 
    (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"]),
    dico_para["tauJ"],dico_para["tauJ"],dico_para["tauM"]]
 end

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
        para_save = (3,false,true) 


# Appelle de la fonction dynamique 
@time module_analyse.heatMapStab(
    model3_ssDD, para_du, para_run, para_matrice_du, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model3_ssDD, para_dl, para_run, para_matrice_dl, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model3_ssDD, para_dm, para_run, para_matrice_dm, para_save, para_plot)
    module_analyse.printTemp()


@time module_analyse.heatMapStab(
    model3_ssDD, para_da, para_run, para_matrice_da, para_save, para_plot)
    module_analyse.printTemp()


