########################################################################
#### Analyse Modèle Parasitoïde Oeufs Sans Compétition (Modèle 1.5) ####
########################################################################

# Charger les packages et les modules 
if Sys.KERNEL == :NT
    include("P:/Article1/fonction/fc_Solveur.jl")
    include("P:/Article1/fonction/fc_dynamique.jl")
else
    include("/Users/Quentin/Desktop/Article1/fonction/fc_dynamique.jl")
    include("/Users/Quentin/Desktop/Article1/fonction/fc_Solveur.jl")
end

using .solveur_maison
using .module_analyse
#using .module_graphique
using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using DataFrames, CSV
using Dates
using Statistics


#########################
## Créaction du modèle ##
#########################
    
# Apport conditionnel pour les premières itérations
function I(t_act)    
    if t_act < (tauU + tauL + tauM) * 2 && t_act > 0
        return 10
    else 
        return 0
    end
end

function bc_model(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - tauL)[7] 
    hist4 = h(p, t - (tauU + tauL + tauM))[4]
    hist5 = h(p, t - (tauL + tauM))[7]
    hist6 = h(p, t - tauJ)[1]
    hist7 = h(p, t - tauJ)[6]
    hist8 = h(p, t - tauU)[6]

    du[1] = (q * u[4] + I(t)) - a * u[1] * u[6] - deltaU * u[1] - (q * hist1 + I(t-tauU)) * exp(u[7])
    du[2] = (q * hist1 + I(t-tauU)) * exp(u[7]) - deltaL * u[2] - (q * hist2 + I(t-tauU-tauL)) * exp(hist3) * exp(- tauL * deltaL)
    du[3] = (q * hist2 + I(t-tauU-tauL)) * exp(hist3) * exp(- tauL * deltaL) - deltaM * u[3] - (q * hist4 + I(t-tauU-tauL-tauM)) * exp(hist5) * exp(- tauL * deltaL) * exp(- tauM * deltaM)
    du[4] = (q * hist4 + I(t-tauU-tauL-tauM)) * exp(hist5) * exp(- tauL * deltaL) * exp(- tauM * deltaM) - deltaA * u[4] 
    du[5] = a * u[1] * u[6] - deltaU * u[5] - a * hist6 * hist7 * exp(- deltaU * tauJ)
    du[6] = a * hist6 * hist7 * exp(- deltaU * tauJ) + I(t) - deltaP * u[6]
    du[7] = - a * (u[6] - hist8) 

end

function bc_model(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p
    
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

function bc_model(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - tauL)[7] 
    hist4 = h(p, t - (tauU + tauL + tauM))[4]
    hist5 = h(p, t - (tauL + tauM))[7]
    hist6 = h(p, t - tauJ)[1]
    hist7 = h(p, t - tauJ)[6]
    hist8 = h(p, t - tauU)[6]

    du[1] = (q * u[4]) - a * u[1] * u[6] - deltaU * u[1] - (q * hist1) * u[7]
    du[2] = (q * hist1) * u[7] - deltaL * u[2] - (q * hist2) * hist3 * exp(- tauL * deltaL)
    du[3] = (q * hist2) * hist3 * exp(- tauL * deltaL) - deltaM * u[3] - (q * hist4 ) * hist5 * exp(- tauL * deltaL) * exp(- tauM * deltaM)
    du[4] = (q * hist4) * hist5 * exp(- tauL * deltaL) * exp(- tauM * deltaM) - deltaA * u[4] 
    du[5] = a * u[1] * u[6] - deltaU * u[5] - a * hist6 * hist7 * exp(- deltaU * tauJ)
    du[6] = a * hist6 * hist7 * exp(- deltaU * tauJ)  - deltaP * u[6]
    du[7] = - a * u[7] * (u[6] - hist8) 

end


####################################################
## Tracé de dynamiques pour le modèle du ravageur ##
####################################################

# Création paramètres et fixation de leur paramètres
#q = 104.6/12.2
#deltaU = 0.0432
q = 70
deltaU = 0.001
deltaL = 0.0048
deltaM = 0.0150
deltaA = 1/12.2
deltaP = 1
tauU = 4.72
tauL = 16
tauM = 8.92
tauJ = tauU
alphaL = 0.01
a = 1

lags = [tauU, (tauU + tauL), tauL, (tauU + tauL + tauM), (tauL + tauM),tauJ,tauJ,tauU]


p = (a, q, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ)


uini = - tauU * (a + deltaU)                # pour Pini = 1
hini = - tauU * deltaU
uStart = 1
u0 = [1, 1, 1, 1, 1, 1, uini]
histini= [1, 1, 1, 1, 1, 1, uini]

tspan = (0.0, 500.0)
alg = "maison"
pdt = 0.001

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
para_run = (u0,histini,tspan,alg,pdt)

# para save 
    # numéro du modèle 
    # présence de compétition chez les larves 
    # sauvgarde ou non 
para_save = (1,false,false) 

# plot 
    # tuple du numéro des variables tracé
    # tuple du nom des variables 
    # taille du segment tracé
para_plot = ((1,2,3,4,5,6),("Oeufs","Larve","Nymphe","Adulte","Oeufs Parasité","Parasitoïde"),500)


# Appelle de la fonction dynamique 
@time gg = module_analyse.traceDynamique(bc_model, para, para_run, para_save, para_plot)
module_analyse.printTemp()


####################################################
## Création de Heatmap pour le modèle de ravageur ## 
####################################################


# Création paramètres et fixation de leur paramètres

# Création des paramètres du modèle
DicoDesParas = Dict(
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
)

tspan = (0.0, 10000.0)
alg = "fixe"
pdt = 0.001


# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("a","q", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")

# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para) 
     return [1, 1, 1, 1, 1, 1, - dico_para["tauU"] * (dico_para["a"] + dico_para["deltaU"])]
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
    return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]), dico_para["tauL"],
    (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauM"]),
    (dico_para["tauL"] + dico_para["tauM"]),dico_para["tauJ"],dico_para["tauJ"],dico_para["tauU"]]
end


# para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
para = (DicoDesParas,fct_lags, ordre_para)

# para_run 
    # fonction créant les valeurs d'initialisation des variables 
    # début et fin de la dynamique 
    # méthode de résolution numérique
para_run = (fct_hini,tspan,alg, pdt)

# para_matrice
    # tuple des deux noms de paramètres qui serviront d'axe à la matrice
    # tuple des valeur min et max des deux paramètres               
    # tuple des dimensions des deux axes de la matrice
    # int des X dernier pas de temps étudié
para_matrice = (("deltaU","q"),(((0.001),(1)),((1),(70))),(5,5),1000)

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (1,false,false) 
    
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
    :"tauJ" => 4.72,
    :"alphaL" =>  0.001,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"deltaU" =>  0.0432,
    :"deltaP" =>  1,
    :"a" => 1,
    :"q" => 106.4/12.2,
)


tspan = (0.0, 10000.0)
alg = "fixe"
pdt = 0.01


# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("a","q", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para) 
    return [1, 1, 1, 1, 1, 1, - dico_para["tauU"] * (dico_para["a"] + dico_para["deltaU"])]
end 

## Fonction perttant de regrouper tout les lags dans un vecteur 
function fct_lags(dico_para)
   return [dico_para["tauU"], (dico_para["tauU"] + dico_para["tauL"]), dico_para["tauL"],
   (dico_para["tauU"] + dico_para["tauL"] + dico_para["tauL"]),
   (dico_para["tauL"] + dico_para["tauL"]),dico_para["tauJ"],dico_para["tauJ"],dico_para["tauU"]]
end

# para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
para = (DicoDesParas,fct_lags, ordre_para)

# para_run 
    # fonction créant les valeurs d'initialisation des variables 
    # début et fin de la dynamique 
    # méthode de résolution numérique
    # pas de temps pour la résolution
para_run = (fct_hini,tspan,alg,pdt)

# para_matrice
    # tuple du nom de paramètres qui serviront d'axe à la matrice
    # tuple des valeur min et max des deux paramètres               
    # nombre de point étudié
    # int des X dernier pas de temps étudié

    para_bifucation = ("deltaA",(0.001,1),3,500)
    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (1,false,false) 
    
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # affiche toutes les valeurs prise ou juste les max/min 
para_plot = (2,true)

# Appelle de la fonction dynamique 

@time rez = module_analyse.simulate_bifurcation(bc_model, para, para_run, para_bifucation, para_save, para_plot)







#######################################################
## Tracer les bifurcations par rapport aux mortalité ##
#######################################################


# Création paramètres et fixation de leur paramètres

# Création des paramètres du modèle

# Création du Dico global
DicoDesParasGlobal = Dict(
    :"tauU" => 4.72,
    :"tauL" => 16,
    :"tauM" => 8.92,
    :"tauJ" => 4.72,
    :"deltaU" =>  0.0432,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"deltaA" => 1/12.2,
    :"deltaP" =>  1,
    :"a" => 1,
    :"q" => 106.4/12.2,
    :"alphaL" => 0.001,
)
# Paramètres de simulation
tspan = (0.0, 10000.0)
alg = "maison"
pdt = 0.005


# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "tauU", "tauL", "tauM", "tau1", "tau2", "tau3", "tau4", "tau5")


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


# para
    # dictionnaire des paramètres du modèle et de leur valeur
    # tuple des valeurs des lags
    # ordre des paramètres dans p 
para = (filter(kv -> kv[1] != "p", DicoDesParasGlobal),fct_lags, ordre_para)  
para_du = (filter(kv -> kv[1] != "deltaU", DicoDesParasGlobal),fct_lags, ordre_para)
para_dl = (filter(kv -> kv[1] != "deltaL", DicoDesParasGlobal),fct_lags, ordre_para)
para_dm = (filter(kv -> kv[1] != "deltaM", DicoDesParasGlobal),fct_lags, ordre_para)
para_da = (filter(kv -> kv[1] != "deltaA", DicoDesParasGlobal),fct_lags, ordre_para)

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
para_bifucation = ("q",(1,70),30,500)
para_bifucation_du = ("deltaU",(0,1),30,500)
para_bifucation_dl = ("deltaL",(0,1),30,500)
para_bifucation_dm = ("deltaM",(0,1),30,500)
para_bifucation_da = ("deltaA",(0,1),30,500)

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (0,false,true) 
    
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # affiche toutes les valeurs prise ou juste les max/min 
para_plot = (2,true)
    
# Appelle de la fonction dynamique 
@time empU = module_analyse.simulate_bifurcation(
    bc_model, para_du, para_run, para_bifucation_du, para_save, para_plot);
if Sys.KERNEL == :NT
    run(`cmd /c start "P:/Article1/graphs/sansPara/temp.png"`)      
else
    run(`cmd /c "/Users/Quentin/Desktop/Article1/graphs/sansPara/temp.png"`)
end 
  
        
@time empL = module_analyse.simulate_bifurcation(
    bc_model, para_dl, para_run, para_bifucation_dl, para_save, para_plot);
if Sys.KERNEL == :NT
    run(`cmd /c start "P:/Article1/graphs/sansPara/temp.png"`)       
else
    run(`cmd /c "/Users/Quentin/Desktop/Article1/graphs/sansPara/temp.png"`)
end 


@time empM = module_analyse.simulate_bifurcation(
    bc_model, para_dm, para_run, para_bifucation_dm, para_save, para_plot);
if Sys.KERNEL == :NT
    run(`cmd /c start "P:/Article1/graphs/sansPara/temp.png"`)       
else
    run(`cmd /c "/Users/Quentin/Desktop/Article1/graphs/sansPara/temp.png"`)
end 


@time empA = module_analyse.simulate_bifurcation(
    bc_model, para_da, para_run, para_bifucation_da, para_save, para_plot);
if Sys.KERNEL == :NT
    run(`cmd /c start "P:/Article1/graphs/sansPara/temp.png"`)       
else
    run(`cmd /c "/Users/Quentin/Desktop/Article1/graphs/sansPara/temp.png"`)
end 





##################################################
## Tracer les heatmap par rapport aux mortalité ##
##################################################


# Création paramètres et fixation de leur paramètres

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


# Paramètres de simulation
tspan = (0.0, 10000.0)
alg = "maison"
pdt = 0.001

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("a", "q", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")

# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        return [0, 0, 0, 0, 0, 0, - dico_para["tauU"] * (dico_para["a"] + dico_para["deltaU"])]
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
para_matrice_du = (("q","deltaU"),((1,70),(0,1)),(2,2),500)
para_matrice_dl = (("q","deltaL"),((1,70),(0,1)),(30,30),500)
para_matrice_dm = (("q","deltaM"),((1,70),(0,1)),(30,30),500)
para_matrice_da = (("q","deltaA"),((1,70),(0,1)),(30,30),500)

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (0,false,true) 
    
        
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # tracer la heatmap de stabilité
        # tracer la heatmap des moyennes 
para_plot = (2,true, true)


# Appelle de la fonction dynamique 
@time module_analyse.heatMapStab(
    bc_model, para_du, para_run, para_matrice_du, para_save, para_plot)
module_analyse.printTemp()


@time module_analyse.heatMapStab(
    bc_model, para_dl, para_run, para_matrice_dl, para_save, para_plot)
module_analyse.printTemp()


@time module_analyse.heatMapStab(
    bc_model, para_dm, para_run, para_matrice_dm, para_save, para_plot)
module_analyse.printTemp()


@time module_analyse.heatMapStab(
    bc_model, para_da, para_run, para_matrice_da, para_save, para_plot)
module_analyse.printTemp()
