#####################################################
#### Analyse Modèle Ravageur 4 stades (Modèle 0) ####
#####################################################

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

    du[1] = q * hist1 * exp(- deltaU * tauU)   -   
        q * hist2 * exp(-deltaU * tauU) * exp(u[3])   -   (deltaL + alphaL * u[1]) * u[1] 
    du[2] = q * hist3 * exp(- deltaU * tauU) * exp(-deltaM * tauM) * exp(hist4)   -   deltaA * u[2] 
    du[3] = - alphaL * (u[1] - hist5) 
end

function model_dev(du, u, h, p, t)
    q, alphaL, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM, tau1, tau2, tau3, tau4, tau5 = p

    hist1 = h(p, t - tau1)[4]
    hist2 = h(p, t - tau2)[4]
    hist3 = h(p, t - tau3)[4] 
    hist4 = h(p, t - tau4)[5]
    hist5 = h(p, t - tau5)[2]

    du[1] = q * u[4] - deltaU * u[1] - q * hist1 * exp(- deltaU * tauU)
    du[2] = q * hist1 * exp(- deltaU * tauU)   -   
        q * hist2 * exp(-deltaU * tauU) * exp(u[5])   -   (deltaL + alphaL * u[2]) * u[2] 
    du[3] = q * hist2 * exp(-deltaU * tauU) * exp(u[5]) - deltaM * u[3] - 
        q * hist3 * exp(- deltaU * tauU) * exp(-deltaM * tauM) * exp(hist4)
    du[4] = q * hist3 * exp(- deltaU * tauU) * exp(-deltaM * tauM) * exp(hist4)   -   deltaA * u[4] 
    du[5] = - alphaL * (u[2] - hist5) 
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
#deltaU = 0
deltaU =  0.0432
deltaL = 0.0048
deltaM = 0.0150
#deltaM = 0.53
#q = 104.6/12.2
q = 70
#deltaA = 1/12.2
deltaA = 0.30
p = (q, alphaL, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM, tau1, tau2, tau3, tau4, tau5)


hini = - tauL * (alphaL + deltaL)
uStart = 1
u0 = [0, uStart, hini]
histini= [0, 0, - tauL * deltaL]


u0_dev = [uStart, uStart, uStart, uStart, hini]
histini_dev= [0, 0 , 0, 0, - tauL * alphaL]

tspan = (0.0, 20000.0)
alg = "fixe"
pdt = 0.01

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
para_run_dev = (u0_dev,histini_dev,tspan,alg,pdt)

# para save 
    # numéro du modèle 
    # présence de compétition chez les larves 
    # sauvgarde ou non 
para_save = (0,true,false) 

# plot 
    # tuple du numéro des variables tracé
    # tuple du nom des variables 
    # taille du segment     racé
    para_plot = ((1),("Larve"),500)


# Appelle de la fonction dynamique 
@time rez = module_analyse.traceDynamique(bc_model, para, para_run, para_save, para_plot);
module_analyse.printTemp()
std(teste[1][end-20000:end])[1]

t_values = rez[2]
u_values = rez[1]

u_obs = u_values[end-500*1000:end]
t_obs = t_values[end-500*1000:end]

# Calculer la variance et la moyenne des 100 derniers pas de temps
mean(u_obs)[1]
std(u_obs)[1]
maximum(u_obs)[1]
minimum(u_obs)[1]

        # Trouver les indices des pics locaux
        peak_indices = findall(i -> u_obs[i-1] < u_obs[i] > u_obs[i+1], 2:length(u_obs)-1) .+ 1  # Corriger les indices

        # Calculer les écarts de temps entre les pics
        if length(peak_indices) > 1
            periods = diff(t_obs[peak_indices])  # Différence entre les instants des pics
            dominant_period = mean(periods)  # Moyenne des périodes détectées
            println(dominant_period)
        
            # Trouver resultsIntegral
            t0 = t_obs[end] - dominant_period  # Prendre la dernière période pour être sûr d’avoir un cycle complet
            indices = findall(x -> t0 ≤ x ≤ t_obs[end], t_obs)
            u1 = [v[1] for v in u_obs]
            integral = trapz(t_obs[indices], u1[indices])
            println(integral)
        end

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

# Paramètres de simulation
tspan = (0.0, 10000.0)

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
para = (DicoDesParas,fct_lags, ordre_para)

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


tspan = (0.0, 20000.0)

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
para = (DicoDesParas,fct_lags, ordre_para)

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
para_bifucation = ("q",(1,70),100,1000)
    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (0,true,true) 
    
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # affiche toutes les valeurs prise ou juste les max/min 
para_plot = (1,true)
    
# Appelle de la fonction dynamique 
@time module_analyse.simulate_bifurcation(bc_model, para, para_run, para_bifucation, para_save, para_plot)






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
#    :"q" => 30,
)

# Paramètres de simulation
tspan = (0.0, 10000.0)
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
para_bifucation = ("q",(1,70),30,50)
para_bifucation_du = ("deltaU",(0,1),30,200)
para_bifucation_dl = ("deltaL",(0,1),30,500)
para_bifucation_dm = ("deltaM",(0,1),30,500)
para_bifucation_da = ("deltaA",(0,1),30,500)

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (0,false,false) 
    
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # affiche toutes les valeurs prise ou juste les max/min 
para_plot = (1,true)
    
# Appelle de la fonction dynamique 
@time empU = module_analyse.simulate_bifurcation(
    bc_model, para_du, para_run, para_bifucation_du, para_save, para_plot);
    module_analyse.printTemp()
  
        
@time empL =module_analyse.simulate_bifurcation(
    bc_model, para_dl, para_run, para_bifucation_dl, para_save, para_plot);
    module_analyse.printTemp()


@time empM = module_analyse.simulate_bifurcation(
    bc_model, para_dm, para_run, para_bifucation_dm, para_save, para_plot);
    module_analyse.printTemp()


@time empA =module_analyse.simulate_bifurcation(
    bc_model, para_da, para_run, para_bifucation_da, para_save, para_plot);
    module_analyse.printTemp() 


#### Création de la fonction de plot bifurcation + equlibre théorique

# Paramètres
#   chm : chemin d'acces au csv
#
#   bif : nom donnée bifurcation
#
#   dfp : (début, fin, pas)

function bifurcationReelle(chm, bif, dfp)

    data = CSV.read("/Users/Quentin/Desktop/Article1/tables/"*chm*".csv", DataFrame; header=false)
    x = data.Column1  
    y = data.Column2  
    plt = plot(x, y, label="Données", xlabel="Colonne 1", ylabel="Colonne 2", title="Bifurcation et valeur d'équilibre")
    axe_range = range(dfp[1],dfp[2],dfp[3])
    for (i, axe) in enumerate(axe_range)
        scatter!(plt,
            fill(axe, length(bif[i])),  # Toutes les valeurs de q (fixe) pour les résultats correspondants
            bif[i],                    # Les premières valeurs extraites
            label=false,                               # Pas de légende
            color=:black,                              # Couleur noire
            alpha=0.4,                                 # Transparence faible (valeurs possibles : 0 à 1)
            ms=2)                                      # Taille des points (marker size)
    end
    display(plt)
    return data
end


bifurcationReelle("sansPara/BifucationU_Mathematica",empU,(para_bifucation_du[2][1],para_bifucation_du[2][2],para_bifucation_du[3]))
bifurcationReelle("sansPara/BifucationL_Mathematica",empL,(para_bifucation_du[2][1],para_bifucation_du[2][2],para_bifucation_du[3]))
bifurcationReelle("sansPara/BifucationM_Mathematica",empM,(para_bifucation_du[2][1],para_bifucation_du[2][2],para_bifucation_du[3]))
bifurcationReelle("sansPara/BifucationA_Mathematica",empA,(para_bifucation_du[2][1],para_bifucation_du[2][2],para_bifucation_du[3]))






##################################################
## Tracer les heatmap par rapport aux mortalité ##
##################################################


# Création paramètres et fixation de leur paramètres

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
tspan = (0.0, 10000.0)
alg = "maison"
pdt = 0.01

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "tauU", "tauL", "tauM", "tau1", "tau2", "tau3", "tau4", "tau5")


# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        return [0, 0, - dico_para["tauL"] * (dico_para["deltaL"])]
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
para_matrice_du = (("q","deltaU"),((1,70),(0.001,1)),(10,10),200)
para_matrice_dl = (("q","deltaL"),((1,70),(0.001,1)),(10,10),200)
para_matrice_dm = (("q","deltaM"),((1,70),(0.001,1)),(10,10),200)
para_matrice_da = (("q","deltaA"),((1,70),(0.001,1)),(10,10),200)

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



# Appelle de la fonction de création des heatmaps 
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



# Création des heatmap en plusieurs couches
module_graphique.heatmap_compose(false, "sansPara", "DD", "q", "deltaU", 200)
module_analyse.printTemp()

module_graphique.heatmap_compose(false, "sansPara", "DD", "q", "deltaL", 200)
module_analyse.printTemp()

module_graphique.heatmap_compose(false, "sansPara", "DD", "q", "deltaM", 200)
module_analyse.printTemp()

module_graphique.heatmap_compose(false, "sansPara", "DD", "q", "deltaA",200, "sansPara/LimiteStabDelta1")
module_analyse.printTemp()




pdt=0.001
# Test resolution, méthod et pdt
@time mean_maison = module_analyse.heatMapStab(
    bc_model, para_du, (fct_hini,tspan,"maison",pdt), para_matrice_du, para_save, para_plot)
module_analyse.printTemp()

@time mean_fixe = module_analyse.heatMapStab(
    bc_model, para_du, (fct_hini,2000,"fixe",pdt), para_matrice_du, para_save, para_plot)
module_analyse.printTemp()

diff = [[mean_maison[1,1] - mean_fixe[1,1],mean_maison[1,2] - mean_fixe[1,2]],[mean_maison[2,1] - mean_fixe[2,1],mean_maison[2,2] - mean_fixe[2,2]]]


