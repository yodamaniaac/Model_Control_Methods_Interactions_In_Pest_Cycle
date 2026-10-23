###############################################
#### Création du diagramme des dynamiques  ####
###############################################


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

model = module_analyse.nymphal_parasitoid_model

####################################################
## Tracé de dynamiques pour le modèle du ravageur ##
####################################################

# Création paramètres et fixation de leur paramètres
DicoDesParasGlobalNul1 = Dict(
    :"deltaU" =>  0,
    :"tauU" => 4.0,
    :"tauL" => 16.0,
    :"tauM" => 8.0,
    :"tauJ" => 16.0,
    :"alphaL" =>  0.0001,
    :"deltaL" => 0.1,
    :"deltaM" => 0.0,
    :"deltaA" => 0.1,
    :"deltaP" =>  0.1,
    :"a" => 0.01,
    :"q" => 30,
    :"H" => 0.25,
    :"mi" => 1,
)

DicoDesParasGlobalNul2 = Dict(
    :"deltaU" =>  0,
    :"tauU" => 4.0,
    :"tauL" => 16.0,
    :"tauM" => 8.0,
    :"tauJ" => 16.0,
    :"alphaL" =>  0.0001,
    :"deltaL" => 0.15,
    :"deltaM" => 0.0,
    :"deltaA" => 0.1,
    :"deltaP" =>  0.1,
    :"a" => 0.01,
    :"q" => 30,
    :"H" => 0.25,
    :"mi" => 1,
)

DicoDesParasGlobalNul3 = Dict(
    :"deltaU" =>  0,
    :"tauU" => 4.0,
    :"tauL" => 16.0,
    :"tauM" => 8.0,
    :"tauJ" => 16.0,
    :"alphaL" =>  0.0001,
    :"deltaL" => 0.2,
    :"deltaM" => 0.0,
    :"deltaA" => 0.1,
    :"deltaP" =>  0.1,
    :"a" => 0.01,
    :"q" => 30,
    :"H" => 0.25,
    :"mi" => 1,
)

DicoDesParasGlobalNul4 = Dict(
    :"deltaU" =>  0,
    :"tauU" => 4.0,
    :"tauL" => 16.0,
    :"tauM" => 8.0,
    :"tauJ" => 16.0,
    :"alphaL" =>  0.0001,
    :"deltaL" => 0.25,
    :"deltaM" => 0.0,
    :"deltaA" => 0.1,
    :"deltaP" =>  0.1,
    :"a" => 0.01,
    :"q" => 30,
    :"H" => 0.25,
    :"mi" => 1,
)

tspan = (0.0, 4000.0)
alg = "fixe"
pdt = 0.0001


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
para1 = (DicoDesParasGlobalNul1,fct_lags, ordre_para)
para2 = (DicoDesParasGlobalNul2,fct_lags, ordre_para)
para3 = (DicoDesParasGlobalNul3,fct_lags, ordre_para)
para4 = (DicoDesParasGlobalNul4,fct_lags, ordre_para)


## Fonction permettant de choisir les valeurs initiales 

function fct_hini(dico_para,type) 
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


# Structurer les paramètres pour les inserer dans le modèle


# para_run 
    # fonction créant les valeurs d'initialisation des variables 
    # début et fin de la dynamique 
    # méthode de résolution numérique
    # int des derniers pas de temps observé

para_run = (fct_hini,tspan,"fixe",500)

# para save 
    # numéro du modèle 
    # Bool save ou non 
    # Save Y or N

para_save = ("Nymphal parasitoid",false,false,3) 


rez = [[],[],[],[]]
# Appelle de la fonction dynamique 
@time rez[1] = module_analyse.runDynamics(model, para1, para_run, para_save);
@time rez[2] = module_analyse.runDynamics(model, para2, para_run, para_save);
@time rez[3] = module_analyse.runDynamics(model, para3, para_run, para_save);
@time rez[4] = module_analyse.runDynamics(model, para4, para_run, para_save);



# Creat Plot 

# para_plot
    # string du titre du plot
    # strig du titre de l'axe des X
    # strig du titre de l'axe des Y
    # int de la dynamique à suivre 
    # vecteur de string du label de chaque courbes

E = [v[1] for v in rez[1][1]]
L = [v[2] for v in rez[1][1]]
N = [v[3] for v in rez[1][1]]
A = [v[4] for v in rez[1][1]]
J = [v[5] for v in rez[1][1]]
P = [v[6] for v in rez[1][1]]
HP = [v[7] for v in rez[1][1]]
HC = [v[8] for v in rez[1][1]]
rez[1][2]

parasitism = N .* P .* DicoDesParasGlobalNul1["a"]
parasitism = parasitism[1:length(parasitism)-200]

enter = DicoDesParasGlobalNul1["q"] .* A[1:length(A)-200] .* exp(-DicoDesParasGlobalNul1["deltaU"] * DicoDesParasGlobalNul1["tauU"]) .* exp.(HC[201:length(HC)])

tx_survie = enter.*exp.(HP[1:length(HP)-200])

Plots.plot(parasitism,label="parasitism")
Plots.plot!(enter,label="enter")
Plots.plot!(N[1:length(N)-200],label="density")
Plots.plot!(exp.(HP[1:length(HP)-200]),label="tx para")
Plots.plot!(tx_survie,label="tx survival")
