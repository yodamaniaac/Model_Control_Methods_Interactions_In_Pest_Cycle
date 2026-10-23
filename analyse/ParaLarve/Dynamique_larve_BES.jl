###############################################
#### Création du diagramme des dynamiqyes  ####
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

model0 = module_analyse.without_parasitoid_model
model1 = module_analyse.eggs_parasitoid_model
model2 = module_analyse.larval_parasitoid_model
model3 = module_analyse.nymphal_parasitoid_model

####################################################
## Tracé de dynamiques pour le modèle du ravageur ##
####################################################

# Création paramètres et fixation de leur paramètres

DicoDesParasGlobalNul_q5_eggs = Dict(
    :"deltaU" =>  0.4,
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
    :"q" => 5,
    :"H" => 0.25,
    :"mi" => 1,
)

DicoDesParasGlobalNul_q40_eggs = Dict(
    :"deltaU" =>  0.4,
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
    :"q" => 40,
    :"H" => 0.25,
    :"mi" => 1,
)

DicoDesParasGlobalNul_q5_larvae = Dict(
    :"deltaU" =>  0.0,
    :"tauU" => 4.0,
    :"tauL" => 16.0,
    :"tauM" => 8.0,
    :"tauJ" => 16.0,
    :"alphaL" =>  0.0001,
    :"deltaL" => 0.125,
    :"deltaM" => 0.0,
    :"deltaA" => 0.1,
    :"deltaP" =>  0.1,
    :"a" => 0.01,
    :"q" => 5,
    :"H" => 0.25,
    :"mi" => 1,
)

DicoDesParasGlobalNul_q40_larvae = Dict(
    :"deltaU" =>  0.0,
    :"tauU" => 4.0,
    :"tauL" => 16.0,
    :"tauM" => 8.0,
    :"tauJ" => 16.0,
    :"alphaL" =>  0.0001,
    :"deltaL" => 0.125,
    :"deltaM" => 0.0,
    :"deltaA" => 0.1,
    :"deltaP" =>  0.1,
    :"a" => 0.01,
    :"q" => 40,
    :"H" => 0.25,
    :"mi" => 1,
)

DicoDesParasGlobalNul_q5_nymphae = Dict(
    :"deltaU" =>  0.0,
    :"tauU" => 4.0,
    :"tauL" => 16.0,
    :"tauM" => 8.0,
    :"tauJ" => 16.0,
    :"alphaL" =>  0.0001,
    :"deltaL" => 0.0,
    :"deltaM" => 0.125,
    :"deltaA" => 0.1,
    :"deltaP" =>  0.1,
    :"a" => 0.01,
    :"q" => 5,
    :"H" => 0.25,
    :"mi" => 1,
)

DicoDesParasGlobalNul_q40_nymphae = Dict(
    :"deltaU" =>  0.0,
    :"tauU" => 4.0,
    :"tauL" => 16.0,
    :"tauM" => 8.0,
    :"tauJ" => 16.0,
    :"alphaL" =>  0.0001,
    :"deltaL" => 0.0,
    :"deltaM" => 0.125,
    :"deltaA" => 0.1,
    :"deltaP" =>  0.1,
    :"a" => 0.01,
    :"q" => 40,
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
para_q5_eggs = (DicoDesParasGlobalNul_q5_eggs ,fct_lags, ordre_para)
para_q40_eggs = (DicoDesParasGlobalNul_q40_eggs ,fct_lags, ordre_para)

para_q5_larvae = (DicoDesParasGlobalNul_q5_larvae ,fct_lags, ordre_para)
para_q40_larvae = (DicoDesParasGlobalNul_q40_larvae ,fct_lags, ordre_para)

para_q5_nymphae = (DicoDesParasGlobalNul_q5_nymphae ,fct_lags, ordre_para)
para_q40_nymphae = (DicoDesParasGlobalNul_q40_nymphae ,fct_lags, ordre_para)

## Fonction permettant de choisir les valeurs initiales 

function fct_hini_0(dico_para,type) 
    if type == "h"
        hStart = 0.0
        hPara = 0.0
        return [hStart, hStart, hStart, hStart,- dico_para["tauL"] * dico_para["alphaL"]]
    else 
        uStart = 1.0
        uPara = 0.1
        return [uStart, uStart, uStart, uStart, - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
    end 
end

function fct_hini_1(dico_para,type) 
    if type == "h"
        hStart = 0.0
        hPara = 0.0
        return [hStart, hStart, hStart, hStart, hPara, hPara, - dico_para["tauU"] * (dico_para["a"] * hPara + dico_para["deltaU"]), 
        - dico_para["tauL"] * (dico_para["alphaL"] * 0 + dico_para["deltaL"])]
    else 
        uStart = 1.0
        uPara = 0.1
        return [0, 0, 0, uStart, 0, uPara, - dico_para["tauU"] * (dico_para["a"] * uPara + dico_para["deltaU"]), 
        - dico_para["tauL"] * (dico_para["alphaL"] * 0 + dico_para["deltaL"])]
    end 
end 

function fct_hini_2(dico_para,type) 
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


# Structurer les paramètres pour les inserer dans le modèle


# para_run 
    # fonction créant les valeurs d'initialisation des variables 
    # début et fin de la dynamique 
    # méthode de résolution numérique
    # int des derniers pas de temps observé
para_run_0 = (fct_hini_0,tspan,"fixe",500)
para_run_1 = (fct_hini_1,tspan,"fixe",500)
para_run_2 = (fct_hini_2,tspan,"fixe",500)
para_run_3 = (fct_hini_3,tspan,"fixe",500)

# para save 
    # numéro du modèle 
    # Bool save ou non 
    # Save Y or N
para_save_0 = ("Without parasitoid",false,false,0)   
para_save_1 = ("Eggs parasitoid",false,false,1) 
para_save_2 = ("Larval parasitoid",false,false,2) 
para_save_3 = ("Nymphal parasitoid",false,false,3) 


# para_plot
    # string du titre du plot
    # strig du titre de l'axe des X
    # strig du titre de l'axe des Y
    # int de la dynamique à suivre 
    # vecteur de string du label de chaque courbes

para_plot = [" ", "Time", "Parasitoid density", 6, ["Without parasitoid"," Eggs parasitoid"," Larval parasitoid","Nymphal parasitoid"]]





### Simulations ###


rez_eggs = [[],[]]
# Appelle de la fonction dynamique 
@time rez_eggs[1] = module_analyse.runDynamics(model1, para_q5_eggs, para_run_1, para_save_1);
@time rez_eggs[2] = module_analyse.runDynamics(model1, para_q40_eggs, para_run_1, para_save_1);

# Creat Plot 
@time graph_eggs = module_analyse.traceDynamics(rez_eggs, para_plot)



rez_larvae = [[],[]]
# Appelle de la fonction dynamique 
@time rez_larvae[1] = module_analyse.runDynamics(model2, para_q5_larvae, para_run_2, para_save_2);
@time rez_larvae[2] = module_analyse.runDynamics(model2, para_q40_larvae, para_run_2, para_save_2);

# Creat Plot 
@time graph_larvae = module_analyse.traceDynamics(rez_larvae, para_plot)



rez_nymphae = [[],[]]
# Appelle de la fonction dynamique 
@time rez_nymphae[1] = module_analyse.runDynamics(model3, para_q5_nymphae, para_run_3, para_save_3);
@time rez_nymphae[2] = module_analyse.runDynamics(model3, para_q40_nymphae, para_run_3, para_save_3);

# Creat Plot 
@time graph_nymphae = module_analyse.traceDynamics(rez_nymphae, para_plot)



## Global Plot 
Plots.plot(graph_eggs, graph_larvae, graph_nymphae; layout = @layout([A; B; C]))