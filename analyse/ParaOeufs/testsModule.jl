#####################################################
#### Analyse Modèle Ravageur 4 stades (Modèle 0) ####
#####################################################

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


##############################################################################
## Créaction du modèle ravageur 4 stades avec competition au stade larvaire ##
##############################################################################

function bc_model(du, u, h, p, t)
    q, alphaL, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM, tau1, tau2, tau3, tau4, tau5 = p

    hist1 = h(p, t - tau1)[2]
    hist2 = h(p, t - tau2)[2]
    hist3 = h(p, t - tau3)[2] 
    hist4 = h(p, t - tau4)[3]
    hist5 = h(p, t - tau5)[1]

    du[1] = q * hist1 * exp(- deltaU * tauU)   -   q * hist2 * exp(-deltaU * tauU) * exp(u[3])   -   (deltaL + alphaL * u[1]) * u[1] 
    du[2] = q * hist3 * exp(- deltaU * tauU) * exp(-deltaM * tauM) * exp(hist4)   -   deltaA * u[2] 
    du[3] = - alphaL * (u[1] - hist5) 
end


####################################################
## Tracé de dynamiques pour le modèle du ravageur ##
####################################################

# Création paramètres et fixation de leur paramètres
tauU = 4.72
tauL = 16
tauM = 8.92
tau1 =  tauU
tau2 = tauU + tauL
tau3 =  tauU + tauL + tauM
tau4 =  tauM
tau5 =  tauL
lags = [tau1, tau2, tau3, tau4, tau5]
alphaL =  0.001
deltaU =  0.0432
deltaL = 0.0048
deltaM = 0.0150
q = 104.6/12.2
deltaA = 1/12.2
p = (q, alphaL, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM, tau1, tau2, tau3, tau4, tau5)


hini = - tauL * (alphaL + deltaL)
uStart = 1
u0 = [uStart, uStart, hini]
histini= [uStart, uStart, hini]

tspan = (0.0, 20000.0)

# Structurer les paramètres pour les inserer dans le modèle

# para
    # tuple des valeurs des paramètres du modèle
    # tuple des valeurs des lags
para = (p,lags)

# para_run 
    # valeur des variable en 0
    # valeur historique des varaibles 
    # début et fin de la dynamique 
    # méthode de résolution numérique
para_run = (u0,histini,tspan,MethodOfSteps(Tsit5()))

# para save 
    # numéro du modèle 
    # présence de compétition chez les larves 
    # sauvgarde ou non 
para_save = (0,true,true) 

# plot 
    # tuple du numéro des variables tracé
    # tuple du nom des variables 
    # taille du segment tracé
para_plot = ((1,2),("Larve","Adulte"),20000)


# Appelle de la fonction dynamique 
module_analyse.traceDynamique(bc_model, para, para_run, para_save, para_plot)










####################################################
## Création de Heatmap pour le modèle de ravageur ## 
####################################################


# Création paramètres et fixation de leur paramètres

# Création des paramètres du modèle
DicoDesParas = Dict(
    :"tauU" => 4.72,
    :"tauL" => 16,
    :"tauM" => 8.92,
    :"tau1" =>  4.72,
    :"tau2" => 4.72 + 16,
    :"tau3" =>  4.72 + 16 + 8.92,
    :"tau4" =>  8.92,
    :"tau5" =>  16,
    :"alphaL" =>  0.001,
    :"deltaU" =>  0.0432,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
)

lags = ["tau1", "tau2", "tau3", "tau4", "tau5"]
tspan = (0.0, 20000.0)

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "tauU", "tauL", "tauM", "tau1", "tau2", "tau3", "tau4", "tau5")

# Création des fonctions 
function fct_hini(dico_para) 
     return [1, 1, - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
end 

# para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
para = (DicoDesParas,lags, ordre_para)

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
    para_matrice = (("deltaA","q"),(((0.001),(1)),((1),(70))),(30,30),1000)

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
    para_save = (0,true,false) 
    
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # tracer la heatmap de stabilité
        # tracer la heatmap des moyennes 
    para_plot = (1,true, true)

# Appelle de la fonction dynamique 
@time module_analyse.heatMapStab(bc_model, para, para_run, para_matrice, para_save, para_plot)













#####################################################################
## Création du diagramme de bifurcation pour le modèle de ravageur ## 
#####################################################################


# Création paramètres et fixation de leur paramètres

# Création des paramètres du modèle
DicoDesParas = Dict(
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
)


lags = ["tau1", "tau2", "tau3", "tau4", "tau5"]
tspan = (0.0, 20000.0)

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "tauU", "tauL", "tauM", "tau1", "tau2", "tau3", "tau4", "tau5")

# Création des fonctions 
function fct_hini(dico_para) 
     return [1, 1, - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
end 

# para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
para = (DicoDesParas,lags, ordre_para)

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
para_bifucation = ("q",(1,70),30,1000)
    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (0,true,false) 
    
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # affiche toutes les valeurs prise ou juste les max/min 
para_plot = (1,false)
    
# Appelle de la fonction dynamique 
@time module_analyse.simulate_bifurcation(bc_model, para, para_run, para_bifucation, para_save, para_plot)