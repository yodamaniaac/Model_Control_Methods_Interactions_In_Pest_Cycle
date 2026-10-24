#####################################################
#### Analyse Modèle Ravageur 4 stades (Modèle 0) ####
#####################################################

# Charger les packages et les modules 


include("fc_dynamique.jl")

using .module_analyse


##############################################################################
## Créaction du modèle ravageur 4 stades avec competition au stade larvaire ##
##############################################################################

function bc_model(du, u, h, p, t)
    q, alphaL, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - tauU - tauL)[4]
    hist3 = h(p, t - tauU - tauL - tauM)[4] 
    hist4 = h(p, t - tauM)[5]
    hist5 = h(p, t - tauL)[2]

    du[1] = q * u[4] - deltaU * u[1] - q * hist1 * exp(- deltaU * tauU) 
    du[2] = q * hist1 * exp(- deltaU * tauU) - (deltaL + alphaL * u[2]) * u[2] - q * hist2 * exp(- deltaU * tauU) * exp(u[5])   
    du[3] = q * hist2 * exp(- deltaU * tauU) * exp(u[5]) - deltaM * u[3] - q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(- deltaM * tauM)
    du[4] = q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(- deltaM * tauM)  -   deltaA * u[4] 
    du[5] = - alphaL * (u[2] - hist5) 
end


##################################################
## Tracer les heatmap par rapport aux mortalité ##
##################################################

# Importer paramètres cluster
typeSimul, pdt, alg, dimHM, tspan, histSave, ptnom = IntegrerARGS(ARGS)

# Création du Dico global
DicoDesParasGlobal = dicoPara("nul")

# Sans Compète
DicoDesParasGlobal["alphaL"] = 0

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "tauU", "tauL", "tauM")


# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        return [0, 0, 0, 0, - dico_para["tauL"] * (dico_para["deltaL"])]
    else 
        return [0, 0, 0, 1, - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
    end 
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
   return [dico_para["tauU"],dico_para["tauU"] + dico_para["tauL"],dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"],dico_para["tauM"],dico_para["tauL"]]
end


# para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
if typeSimul == 1
    para_du = (filter(kv -> kv[1] !="deltaU", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)
elseif typeSimul == 2
    para_dl = (filter(kv -> kv[1] !="deltaL", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)
elseif typeSimul == 3
    para_dm = (filter(kv -> kv[1] !="deltaM", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)
elseif typeSimul == 4
    para_da = (filter(kv -> kv[1] !="deltaA", filter(kv -> kv[1] !="q", DicoDesParasGlobal)),fct_lags, ordre_para)
end

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
    
dimPara = ((1,70),(0.001,1))    
dimPara_da = ((1,70),(0.1,1.1))  
dimHeatMap = (dimHM,dimHM)
tstud = 500

if typeSimul == 1
    para_matrice_du = (("q","deltaU"),dimPara,dimHeatMap,tstud)
elseif typeSimul == 2
    para_matrice_dl = (("q","deltaL"),dimPara,dimHeatMap,tstud)
elseif typeSimul == 3
    para_matrice_dm = (("q","deltaM"),dimPara,dimHeatMap,tstud)
elseif typeSimul == 4
    para_matrice_da = (("q","deltaA"),dimPara_da,dimHeatMap,tstud)
end


    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (0,true,true, histSave,ptnom) 
    
        
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # tracer la heatmap de stabilité
        # tracer la heatmap des moyennes 
para_plot = ([1,2,3,4],true, true)


# Appelle de la fonction de création des heatmaps 

if typeSimul == 1
    println("Fait")
    module_analyse.heatMapStab(bc_model, para_du, para_run, para_matrice_du, para_save, para_plot)
elseif typeSimul == 2
    println("Fait")
    module_analyse.heatMapStab(bc_model, para_dl, para_run, para_matrice_dl, para_save, para_plot)
elseif typeSimul == 3
    println("Fait")
    module_analyse.heatMapStab(bc_model, para_dm, para_run, para_matrice_dm, para_save, para_plot)
elseif typeSimul == 4
    println("Fait")
    module_analyse.heatMapStab(bc_model, para_da, para_run, para_matrice_da, para_save, para_plot)
end


