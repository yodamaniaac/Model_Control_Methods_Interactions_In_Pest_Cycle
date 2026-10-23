########################################
#### Diagramme de bifurcation Final ####
########################################

# Charger les packages et les modules 
include("/Users/Quentin/Desktop/Article1/fonction/fc_dynamique.jl")
include("/Users/Quentin/Desktop/Article1/fonction/fc_Solveur.jl")

using .solveur_maison
using .module_analyse
using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using DataFrames, CSV
using Dates

using Statistics
using ColorSchemes
using CairoMakie

#########################
## Créaction du modèle ##
#########################


model3 = module_analyse.nymphal_parasitoid_model


#################################################
## Création du diagramme de bifurcation Finaux ## 
#################################################

#### Paramètres Globaux ####

# Création paramètres et fixation de leur paramètres

# Création des paramètres du modèle
DicoDesParasGlobalNul = Dict(
    :"deltaU" =>  0,
    :"tauU" => 4.0,
    :"tauL" => 16.0,
    :"tauM" => 8.0,
    :"tauJ" => 16.0,
    :"alphaL" =>  0.0001,
    :"deltaL" => 0.0,
    :"deltaM" => 0.0,
    :"deltaA" => 0.1,
    :"deltaP" =>  0.1,
    :"a" => 0.01,
    :"q" => 30,
    :"H" => 0.25,
    :"mi" => 1,
)

tspan = (0.0, 30000.0)

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "a", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ","H", "mi")

## Fonction pertant de regrouper tout les lags dans un vecteur 
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
para = (DicoDesParasGlobalNul,fct_lags, ordre_para)

# valeur maximum de larve sans parasitoïde
#max = 4942.0 # Pour q=70
#max = 3411.45 # Pour q=20
max = 1.0 

# para_matrice
    # tuple des deux noms de paramètres qui serviront d'axe à la matrice
    # tuple des valeur min et max des deux paramètres               
    # tuple des dimensions des deux axes de la matrice
    # int des X dernier pas de temps étudié
para_bifucation = (["deltaL"],(0.00001,0.4),100,5000,([" "]),"ratio")



#### Paramètre Spécifique ####

## Fonction permettant de choisir les valeurs initiales 

function fct_hini_3(dico_para,type) 
    if type == "h"
        hStart = 0
        hPara = 0

        return [hStart, hStart, hStart, hStart, hPara, hPara, 
        - dico_para["tauM"] * (dico_para["a"] * hPara + dico_para["deltaM"]),
        - dico_para["tauL"] * (dico_para["alphaL"] * 0 + dico_para["deltaL"])]
    else 
        uStart = 1
        uPara = 0.1
        return [0, 0, 0, uStart, 0, uPara, 
        - dico_para["tauM"] * (dico_para["a"] * uPara + dico_para["deltaM"]),
        - dico_para["tauL"] * (dico_para["alphaL"] * 0 + dico_para["deltaL"])]
    end 
end 



# para_run 
    # fonction créant les valeurs d'initialisation des variables 
    # début et fin de la dynamique 
    # méthode de résolution numérique
para_run_3 = (fct_hini_3,tspan)


# para save 
    # numéro du modèle 
    # Bool save ou non 
    # Save Y or N
para_save_3 = ("Nymphal parasitoid",false,false,3) 

# type 
    # Type de bifurcation 
type = ["bifurcation_mortalite"]


    

popObs = 7
@time rezRatio = module_analyse.simulate_bifurcation_Final(model3, para, para_run_3, para_bifucation, para_save_3, max, type, popObs);

# Créer les plots 
para_plot_panel = "simple"


# para_plot
    # vecteur de float, limité minimum et maximum de l'axe Y
    # string du nom de l'axe X
    # string du nom de l'axe
para_plot = [[0.0,100.0],"Competition intensity","Survival rate of nymphs"]


module_analyse.plot_bifurcation_Final(rezRatio, para_bifucation, para_save_3, para_plot_panel, para_plot)
