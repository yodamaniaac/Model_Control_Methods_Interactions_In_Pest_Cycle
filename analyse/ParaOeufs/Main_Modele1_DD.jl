######################################################################
#### Analyse Modèle Parasitoïde Oeufs avec Compétition (Modèle 1) ####
######################################################################

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

function bc_model(du, u, h, p, t)
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


####################################################
## Tracé de dynamiques pour le modèle du ravageur ##
####################################################

# Création paramètres et fixation de leur paramètres
#q = 104.6/12.2
deltaU =  0.1
q = 20 
#deltaL = 0.001
deltaL = 0.1
deltaM = 0.1
deltaA = 0.1
deltaP = 0.1
tauU = 4
tauL = 16
tauM = 8
tauJ = tauU
a = 0.001
alphaL = 0.0001


tspan = (0.0, 10000.0)
alg = "fixe"
pdt = 0.0001

lags = [tauU, (tauU + tauL), tauL, (tauU + tauL + tauM), (tauL + tauM), tauM , tauJ, tauJ, tauU, tauL]

p = (q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ)


uStart = 1
uPara = 0.1
hStart = 0
hPara = 0


uini1 = - tauU * (uPara * a + deltaU)
uini2 = - tauL * (uStart * alphaL + deltaL)

hini1 = - tauU * (hPara * a + deltaU)
hini2 = - tauL * (hStart * alphaL + deltaL)

u0 = [0, 0, 0, uStart, 0, uPara, uini1, uini2]
histini= [0, 0, 0, hStart, 0, hPara, hini1, hini2]


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
para_run = (u0,histini,tspan,alg, pdt)
# para save 
    # numéro du modèle 
    # présence de compétition chez les larves 
    # sauvgarde ou non 

para_save = (1,true,true) 

# plot 
    # tuple du numéro des variables tracé
    # tuple du nom des variables 
    # taille du segment tracé
para_plot = ((1,2,3,4,5,6),("Eggs","Larva","Pupa","Adult","Parasitised Eggs","Adult Parasitoid"),500)


# Appelle de la fonction dynamique 
@time rez = module_analyse.traceDynamique(bc_model, para, para_run, para_save, para_plot);
module_analyse.printTemp()

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
    :"tauJ" => 4.72,
    :"alphaL" =>  0.001,
    :"deltaL" => 0.0048,
    :"deltaM" => 0.0150,
    :"deltaA" => 1/12.2,
    :"deltaP" =>  1,
    :"a" => 1,
    :"alphaL" => 0.001,
)

tspan = (0.0, 20000.0)

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "a", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")


# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        hStart = 0
        hPara = 0
        return [hStart, hStart, hStart, hStart, hPara, hPara, - dico_para["tauU"] * 
        (dico_para["a"] * hPara + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] * hStart + dico_para["deltaL"])]
    else 
        uStart = 1
        uPara = 0.1
        return [uStart, uStart, uStart, uStart, uPara, uPara, - dico_para["tauU"] * 
        (dico_para["a"] * uPara + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] * uStart + dico_para["deltaL"])]
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
para_matrice = (("deltaU","q"),(((0.001),(1)),((1),(70))),(5,5),1000)

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (1,true,false) 
    
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # tracer la heatmap de stabilité
        # tracer la heatmap des moyennes 
para_plot = (1,true, true)

# Appelle de la fonction dynamique 
@time module_analyse.heatMapStab(bc_model, para, para_run, para_matrice, para_save, para_plot)
module_analyse.printTemp()


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
    :"deltaA" => 1/12.2,
    :"deltaP" =>  1,
    :"a" => 1,
    :"q" => 106.4/12.2,
    :"alphaL" => 0.001,
)


tspan = (0.0, 20000.0)

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")


# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para) 
    return [1, 1, 1, 1, 1, 1, - dico_para["tauU"] * (dico_para["a"] + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] + dico_para["deltaL"])]
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
para_bifucation = ("deltaU",(0.001,1),30,1000)
    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (1,true,false) 
    
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
    :"tauJ" => 4.72,
    :"alphaL" =>  0.001,
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
alg = "fixe"
pdt = 0.001


# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q","a", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")


# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        hStart = 0
        hPara = 0
        return [hStart, hStart, hStart, hStart, hPara, hPara, - dico_para["tauU"] * 
        (dico_para["a"] * hPara + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] * hStart + dico_para["deltaL"])]
    else 
        uStart = 1
        uPara = 0.1
        return [uStart, uStart, uStart, uStart, uPara, uPara, - dico_para["tauU"] * 
        (dico_para["a"] * uPara + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] * uStart + dico_para["deltaL"])]
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
para_bifucation = ("q",(1,70),5,200)
para_bifucation_du = ("deltaU",(0,1),5,200)
para_bifucation_dl = ("deltaL",(0,1),30,500)
para_bifucation_dm = ("deltaM",(0,1),30,500)
para_bifucation_da = ("deltaA",(0,1),30,500)

    # para save 
        # numéro du modèle 
        # Bool présence de densité dépendence sur les larves 
        # Bool save ou non 
para_save = (0,true,true) 
    
    # plot 
        # tuple du numéro de la variables dont la stabilité est observé
        # affiche toutes les valeurs prise ou juste les max/min 
para_plot = (2,true)
    
# Appelle de la fonction dynamique 
@time empU = module_analyse.simulate_bifurcation(
    bc_model, para_du, para_run, para_bifucation_du, para_save, para_plot);
module_analyse.printTemp()
  
        
@time empL = module_analyse.simulate_bifurcation(
    bc_model, para_dl, para_run, para_bifucation_dl, para_save, para_plot);
module_analyse.printTemp()

@time empM = module_analyse.simulate_bifurcation(
    bc_model, para_dm, para_run, para_bifucation_dm, para_save, para_plot);
module_analyse.printTemp()

@time empA = module_analyse.simulate_bifurcation(
    bc_model, para_da, para_run, para_bifucation_da, para_save, para_plot);
module_analyse.printTemp()





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
)


# Paramètres de simulation
tspan = (0.0, 20000.0)
alg = "fixe"
pdt = 0.0001

# Structurer les paramètres pour les inserer dans le modèle
ordre_para = ("q","a", "alphaL", "deltaU", "deltaL", "deltaM", "deltaA", "deltaP", "tauU", "tauL", "tauM", "tauJ")


# Création des fonctions 

## Fonction permettant de choisir les valeurs initiales 
function fct_hini(dico_para,type) 
    if type == "h"
        hStart = 0
        hPara = 0
        return [hStart, hStart, hStart, hStart, hPara, hPara, - dico_para["tauU"] * 
        (dico_para["a"] * hPara + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] * hStart + dico_para["deltaL"])]
    else 
        uStart = 1
        uPara = 0.1
        return [0, 0, 0, uStart, 0, uPara, - dico_para["tauU"] * 
        (dico_para["a"] * uPara + dico_para["deltaU"]),
        - dico_para["tauL"] * (dico_para["alphaL"] * uStart + dico_para["deltaL"])]
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
dimHeatMap = (5,5)
tstud = 200

para_matrice_du = (("q","deltaU"),dimPara,dimHeatMap,tstud)
para_matrice_dl = (("q","deltaL"),dimPara,dimHeatMap,tstud)
para_matrice_dm = (("q","deltaM"),dimPara,dimHeatMap,tstud)
para_matrice_da = (("q","deltaA"),dimPara,dimHeatMap,tstud)

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









pdt=0.001
# Test resolution, méthod et pdt
@time mean_maison = module_analyse.heatMapStab(
    bc_model, para_du, (fct_hini,tspan,"maison",pdt), para_matrice_du, para_save, para_plot)
#module_analyse.printTemp()

@time mean_fixe = module_analyse.heatMapStab(
    bc_model, para_du, (fct_hini,2000,"fixe",pdt), para_matrice_du, para_save, para_plot)
#module_analyse.printTemp()

diff = [[mean_maison[1,1] - mean_fixe[1,1],mean_maison[1,2] - mean_fixe[1,2]],[mean_maison[2,1] - mean_fixe[2,1],mean_maison[2,2] - mean_fixe[2,2]]]
