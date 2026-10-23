#######################################################################
#### Analyse Modèle Parasitoïde larves avec Compétition (Modèle 2) ####
#######################################################################

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
    hist3 = h(p, t - (tauU + tauL + tauM))[4]
    hist4 = h(p, t - tauM)[7]
    hist5 = h(p, t - tauJ)[2]
    hist6 = h(p, t - tauJ)[6]
    hist7 = h(p, t - tauL)[6]
    hist8 = h(p, t - tauL)[2]
    hist9 = h(p, t - tauL)[5]
    hist10 = h(p, t - tauJ)[2]
    hist11 = h(p, t - tauJ)[5]


    du[1] = q * u[4] - deltaU * u[1] - q * hist1 * exp(- deltaU * tauU)
    du[2] = q * hist1 * exp(- deltaU * tauU) - a * u[2] * u[6] - 
        (alphaL * (u[2] + u[5]) + deltaL) * u[2] - q * hist2 * exp(- deltaU * tauU) * exp(u[7])
    du[3] = q * hist2 * exp(- deltaU * tauU) * exp(u[7]) - deltaM * u[3] - 
    q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(- tauM * deltaM)
    du[4] = q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(- tauM * deltaM)   -   deltaA * u[4] 
    du[5] = a * u[2] * u[6] - (alphaL * (u[2] + u[5]) + deltaL) * u[5] - a * hist5 * hist6 * exp(u[8])
    du[6] = a * hist5 * hist6 * exp(u[8]) - deltaP * u[6]
    du[7] = - a * (u[6] - hist7) - alphaL * ((u[2] + u[5]) - (hist8 + hist9))
    du[8] = - alphaL * ((u[2] + u[5])-(hist10 + hist11))
end


####################################################
## Tracé de dynamiques pour le modèle du ravageur ##
####################################################

# Création paramètres et fixation de leur paramètres
#q = 104.6/12.2
deltaU =  0
q = 70 
#deltaL = 0.001
deltaL = 0.0
deltaM = 0.0
deltaA = 0.1
deltaP = 0.1
tauU = 4
tauL = 16
tauM = 8
tauJ = tauL
a = 0.01
alphaL = 0.00001


tspan = (0.0, 10000.0)
alg = "fixe"
pdt = 0.0001

lags = [tauU, (tauU + tauL), tauL, (tauU + tauL + tauM), (tauL + tauM), tauM , tauJ, tauJ, tauU, tauL]

p = (q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ)


uStart = 1
uPara = 0.1
hStart = 0
hPara = 0


uini1 = - tauL * (uPara * alphaL + deltaL + a * uPara)
uini2 = - tauJ * (uPara * alphaL + deltaL)

hini1 = - tauL * (hPara * alphaL + deltaL + a * hPara)
hini2 = - tauJ * (hPara * alphaL + deltaL)

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
@time rez2 = module_analyse.traceDynamique(bc_model, para, para_run, para_save, para_plot);
module_analyse.printTemp()



mean([x[1] for x in rez2[1][end-1000000:end]])
mean([x[2] for x in rez2[1][end-1000000:end]])
mean([x[3] for x in rez2[1][end-1000000:end]])
mean([x[4] for x in rez2[1][end-1000000:end]])
mean([x[5] for x in rez2[1][end-1000000:end]])
mean([x[6] for x in rez2[1][end-1000000:end]])
