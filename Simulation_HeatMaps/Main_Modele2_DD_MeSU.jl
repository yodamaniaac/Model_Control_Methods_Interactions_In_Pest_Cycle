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
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ, H = p

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

function bc_model_II(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ, H = p

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
    du[2] = q * hist1 * exp(- deltaU * tauU) - (a * u[2] * u[6]) / (1 + u[2] * a * H) - 
       (alphaL * (u[2] + u[5]) + deltaL) * u[2] - q * hist2 * exp(- deltaU * tauU) * exp(u[7])
    du[3] = q * hist2 * exp(- deltaU * tauU) * exp(u[7]) - deltaM * u[3] - 
    q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(- tauM * deltaM)
    du[4] = q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(- tauM * deltaM)   -   deltaA * u[4] 
    du[5] = (a * u[2] * u[6]) / (1 + u[2] * a * H) - (alphaL * (u[2] + u[5]) + deltaL) * u[5] - (a * hist5 * hist6) / (1 + hist5 * a * H) * exp(u[8])
    du[6] = (a * hist5 * hist6) / (1 + hist5 * a * H) * exp(u[8]) - deltaP * u[6]
    du[7] = - (a * (u[6] - hist7)) / (1 + a * H * (u[2] - hist8)) - alphaL * ((u[2] + u[5]) - (hist8 + hist9))
    du[8] = - alphaL * ((u[2] + u[5])-(hist10 + hist11))
end



##################################################
## Tracer les heatmap par rapport aux mortalité ##
##################################################

# Importer paramètres cluster
typeSimul, pdt, alg, dimHM, tspan, histSave, ptnom = IntegrerARGS(ARGS)


# Création paramètres et fixation de leur paramètres

# Création des paramètres du modèle

# Création du Dico global
DicoDesParasGlobal = dicoPara("nul")


# Si dico des para normal 
#DicoDesParasGlobal["tauJ"] = DicoDesParasGlobal["tauL"]
DicoDesParasGlobal["tauJ"] = 16

# Si RF type 2, h différentde 0
DicoDesParasGlobal["H"] = 0



# Si dico des para normal 
#DicoDesParasGlobal["tauJ"] = 16

# Sans Compète
DicoDesParasGlobal["alphaL"] = 0

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "a", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ", "H")



# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        hStart = 0
        hPara = 0
        return [hStart, hStart, hStart, hStart, hPara, hPara, 
        - dico_para["tauL"] * ((hPara + 0) * dico_para["alphaL"] + dico_para["deltaL"] + 
        dico_para["a"] * hPara),
        - dico_para["tauJ"] * ((hPara + 0) * dico_para["alphaL"] + dico_para["deltaL"])]
    else 
        uPara = 0.1
        uStart = 1
        return [0, 0, 0, uStart, 0, uPara, 
        - dico_para["tauL"] * ((uPara + 0) * dico_para["alphaL"] + dico_para["deltaL"] + 
        dico_para["a"]* uPara),
        - dico_para["tauJ"] * ((uPara + 0) * dico_para["alphaL"] + dico_para["deltaL"])]
    end 
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
tstud = 1000

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
para_save = (2,true,true, histSave, ptnom)  
    
        
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # tracer la heatmap de stabilité
        # tracer la heatmap des moyennes 
para_plot = ([1,2,3,4,5,6],true, true)


# Appelle de la fonction de création des heatmaps 

if typeSimul == 1
    println("Fait")
    module_analyse.heatMapStab(bc_model_II, para_du, para_run, para_matrice_du, para_save, para_plot)
elseif typeSimul == 2
    println("Fait")
    module_analyse.heatMapStab(bc_model_II, para_dl, para_run, para_matrice_dl, para_save, para_plot)
elseif typeSimul == 3
    println("Fait")
    module_analyse.heatMapStab(bc_model_II, para_dm, para_run, para_matrice_dm, para_save, para_plot)
elseif typeSimul == 4
    println("Fait")
    module_analyse.heatMapStab(bc_model_II, para_da, para_run, para_matrice_da, para_save, para_plot)
end





