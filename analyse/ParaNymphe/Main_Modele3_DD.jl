######################################################################
#### Analyse Modèle Parasitoïde larves avec Compétition (Modèle 2) ###
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

function model(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - (tauU + tauL + tauM))[4]
    hist4 = h(p, t - tauM)[8]
    hist5 = h(p, t - tauJ)[3]
    hist6 = h(p, t - tauJ)[6]
    hist7 = h(p, t - tauM)[6]
    hist8 = h(p, t - tauL)[2]


    du[1] = q * u[4] - deltaU * u[1] - q * hist1 * exp(- deltaU * tauU) 
    du[2] = q * hist1 * exp(- deltaU * tauU) - (alphaL * u[2] + deltaL) * u[2] -
        q * hist2 * exp(- deltaU * tauU) * exp(u[8])
    du[3] = q * hist2 * exp(- deltaU * tauU) * exp(u[8]) - a * u[3] * u[6] - deltaM * u[3] - 
        q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(u[7])
    du[4] = q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(u[7])   -   deltaA * u[4] 
    du[5] = a * u[3] * u[6] - deltaM * u[5] - a * hist5 * hist6 * exp(- deltaM * tauJ)
    du[6] = a * hist5 * hist6 * exp(- deltaM * tauJ) - deltaP * u[6]
    du[7] = - a * (u[6] - hist7) 
    du[8] = - alphaL * (u[2] - hist8)
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


uini1 = - tauM * (a * uPara + deltaM)
uini2 = - tauL * (alphaL * 0 + deltaL)

hini1 = - tauM * (a * hPara + deltaM)
hini2 = - tauL * (alphaL * 0 + deltaL)

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