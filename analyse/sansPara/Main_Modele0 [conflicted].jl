#####################################################
#### Analyse Modèle Ravageur 4 stades (Modèle 0) ####
#####################################################

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
#deltaM = 0.0150
deltaM = 0.53
q = 104.6/12.2
#q = 30
deltaA = 1/12.2
#deltaA = 0.75
p = (q, alphaL, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM, tau1, tau2, tau3, tau4, tau5)


hini = - tauL * (alphaL + deltaL)
uStart = 1
u0 = [uStart, uStart, hini]
histini= [0, 0, - tauL * alphaL]
u0_dev = [uStart, uStart, uStart, uStart, hini]
histini_dev= [0, 0 , 0, 0, - tauL * alphaL]

tspan = (0.0, 10000.0)
alg = "maison"
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
    # taille du segment tracé
para_plot = ((1,2),("Larve","Adulte"),500)


# Appelle de la fonction dynamique 
@time teste = module_analyse.traceDynamique(bc_model, para, para_run, para_save, para_plot);

plot(teste[1],teste[2])

    
    # Extraire les paramètres de sauvgarde
    nb_modele = para_save[1]
    presence_dd = para_save[2]
    YN_save = para_save[3]

    # Extraire les paramètres pour le graphique
    pop_plot = para_plot[1]
    nom_pop_plot = para_plot[2]
    t_plot = para_plot[3]

    # Extraire les paramètres du modèle 
    p = para[1]
    lags =  para[2]

    # Extraire les paramètres de simulation
    u0 = para_run[1]
    h(p,t) = para_run[2]
    h0 = para_run[2]
    tspan = para_run[3]
    alg = para_run[4]
    if length(para_run) >= 5
        pdt = para_run[5]
    else 
        pdt = 0.01
    end

    # Fixer l'horraire
    heure = string(now())

    modele = bc_model

    # Résoudre le problème selon le solveur choisis 
        # Résoudre le problème selon le solveur choisis 
        if alg == "maison"
            # Définir le problème et le résoudre
            t_values, u_values = solveur_maison.solve_dde_rk2(modele, lags, tspan, u0, pdt, p, h0)
        elseif alg == "fixe"
              # Définir le problème et le résoudre
              prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
              sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=pdt)  
       
              # Extraire les solutions
              t_values = sol.t # Les temps où la solution est évaluée
              u_values = sol.u # Les valeurs des solutions correspondantes
        elseif alg == "hard"
            # Définir le problème et le résoudre
            prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
            sol = solve(prob, MethodOfSteps(Vern9()), abstol=1e-9, reltol=1e-9, dt=1e-5)

            #  Extraire les solutions
            t_values = sol.t # Les temps où la solution est évaluée
            u_values = sol.u # Les valeurs des solutions correspondantes
        else 
            # Définir le problème et le résoudre
            prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
            sol = solve(prob, MethodOfSteps(Tsit5()))  
     
            # Extraire les solutions
            t_values = sol.t # Les temps où la solution est évaluée
            u_values = sol.u # Les valeurs des solutions correspondantes
        end
        t_values
        u_values

    # Séléctionner les valeurs tracé
    last_values_indices = max(length(t_values) - t_plot, 1):length(t_values)  # Assure qu'on ne dépasse pas les limites
    t_last_values = [x for x in t_values if x > (tspan[2] - t_plot)]
    u_last_values = u_values[end-length(t_last_values)+1:end] 

    t_lite = t_last_values[1:Int(round((1/pdt))):end]
    u_lite = u_last_values[1:Int(round((1/pdt))):end]

    # Créer le plot vide
    plt = plot(title="Dynamique du système", xlabel="Temps", ylabel="Densité des stades", lw=2)

    # Ajouter les courbes
    for i in pop_plot
        plot!(plt, t_lite, [u[i] for u in u_lite], label=nom_pop_plot[i])
    end

    # Afficher le tracé
    display(plt)   

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
para_bifucation_du = ("deltaU",(0,1),30,500)
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
    bc_model, para_du, para_run, para_bifucation_du, para_save, para_plot)
@time empL =module_analyse.simulate_bifurcation(
    bc_model, para_dl, para_run, para_bifucation_dl, para_save, para_plot)

@time empM = module_analyse.simulate_bifurcation(
    bc_model, para_dm, para_run, para_bifucation_dm, para_save, para_plot)

@time empA =module_analyse.simulate_bifurcation(
    bc_model, para_da, para_run, para_bifucation_da, para_save, para_plot)

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
para_matrice_du = (("q","deltaU"),((1,70),(0,1)),(5,5),1000)
para_matrice_dl = (("q","deltaL"),((1,70),(0,1)),(5,5),1000)
para_matrice_dm = (("q","deltaM"),((1,70),(0,1)),(5,5),1000)
para_matrice_da = (("q","deltaA"),((1,70),(0,1)),(5,5),1000)

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
@time module_analyse.heatMapStab(
    bc_model, para_du, para_run, para_matrice_du, para_save, para_plot)
@time module_analyse.heatMapStab(
    bc_model, para_dl, para_run, para_matrice_dl, para_save, para_plot)
@time module_analyse.heatMapStab(
    bc_model, para_dm, para_run, para_matrice_dm, para_save, para_plot)
@time module_analyse.heatMapStab(
    bc_model, para_da, para_run, para_matrice_da, para_save, para_plot)


