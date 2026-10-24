
##############################################
#### Diagramme survie du parasitoid en 1D ####
##############################################

# Charger les packages et les modules
include(joinpath(@__DIR__, "fc_local_simulation.jl"))

using .module_analyse
using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using DataFrames, CSV, DelimitedFiles
using Dates

using Statistics
using ColorSchemes
using CairoMakie

###########################
## Créaction des modèles ##
###########################

model2 = module_analyse.larval_parasitoid_model


#################################################
## Création du diagramme de bifurcation Finaux ## 
#################################################

#### Paramètres Globaux ####

# Création paramètres et fixation de leur paramètres

# Création des paramètres du modèle
DicoDesParasGlobalNul_q70 = Dict(
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
    :"q" => 70,
    :"H" => 0.25,
    :"mi" => 1,
)


DicoDesParasGlobalNul_q20 = Dict(
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
    :"q" => 20,
    :"H" => 0.25,
    :"mi" => 1,
)

tspan = (0.0, 10000.0)

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
para_q20 = (DicoDesParasGlobalNul_q20,fct_lags, ordre_para)
para_q70 = (DicoDesParasGlobalNul_q70,fct_lags, ordre_para)

# valeur maximum de larve sans parasitoïde
#max = 4942.0 # Pour q=70
#max = 3411.45 # Pour q=20
max = 1.0 

# para_matrice
    # tuple des deux noms de paramètres qui serviront d'axe à la matrice
    # tuple des valeur min et max des deux paramètres               
    # tuple des dimensions des deux axes de la matrice
    # int des X dernier pas de temps étudié
para_bifucation = ("alphaL",(0.00001,0.001),100,1000,
    ("Without Mortality","Egg Mortality","Larvae Mortality", "Nymphae Mortality"),"log")



    #### Paramètre Spécifique ####

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

# para_run 
    # fonction créant les valeurs d'initialisation des variables 
    # début et fin de la dynamique 
    # méthode de résolution numérique
para_run = (fct_hini,tspan)

# para save 
    # numéro du modèle 
    # Bool save ou non 
    # Save Y or N
para_save = ("Larval parasitoid",false,true,2) 

# type 
    # Type de bifurcation 
type = ["1D"]

#### Run les simulations ####

@time rez1D_q70 = module_analyse.simulate_bifurcation_Final(model2, para_q70, para_run, para_bifucation, para_save, max, type);
@time rez1D_q20 = module_analyse.simulate_bifurcation_Final(model2, para_q20, para_run, para_bifucation, para_save, max, type);
module_analyse.Graph_survie_para_1D([rez1D_q70, rez1D_q20],para_bifucation,false,false, true)
