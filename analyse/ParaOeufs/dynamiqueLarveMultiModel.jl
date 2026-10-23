#########################################################################
#### Analyse Dynamique Différents Modèle Parasitoïde avec Compétition ###
#########################################################################


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
using ColorSchemes


###########################
## Créaction des modèles ##
###########################

function bc_model0(du, u, h, p, t)
    q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU+tauL))[4]
    hist3 = h(p, t - (tauU + tauL + tauM))[4] 
    hist4 = h(p, t - tauM)[5]
    hist5 = h(p, t - tauL)[2]

    du[1] = q * u[4] - deltaU * u[1] - q * hist1 * exp(- deltaU * tauU)
    du[2] = q * hist1 * exp(- deltaU * tauU)   -   
        q * hist2 * exp(-deltaU * tauU) * exp(u[5])   -   (deltaL + alphaL * u[2]) * u[2] 
    du[3] = q * hist2 * exp(-deltaU * tauU) * exp(u[5]) - deltaM * u[3] - 
        q * hist3 * exp(- deltaU * tauU) * exp(-deltaM * tauM) * exp(hist4)
    du[4] = q * hist3 * exp(- deltaU * tauU) * exp(-deltaM * tauM) * exp(hist4)   -   deltaA * u[4] 
    du[5] = - alphaL * (u[2] - hist5) 
end

function bc_model1(du, u, h, p, t)
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

function bc_model2(du, u, h, p, t)
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

function bc_model3(du, u, h, p, t)
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
deltaU =  0
q = 20
deltaL = 0
deltaM = 0
deltaA = 0.1
deltaP = 0.1
tauU = 4
tauL = 16
tauM = 8
tauJ = 16
a = 0.01
alphaL = 0.0001


tspan = (0.0, 3000.0)
alg = "fixe"
pdt = 0.0001

lags = [tauU, (tauU + tauL), tauL, (tauU + tauL + tauM), (tauL + tauM), tauM , tauJ, tauJ, tauU, tauL]

p = (q, a, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ)


uStart = 1
uPara = 0.1
hStart = 0
hPara = 0


u00 = [uStart, uStart, uStart, uStart, - tauL * (alphaL + deltaL)]
histini0= [0, 0 , 0, 0, - tauL * alphaL]

u01 = [0, 0, 0, uStart, 
    0, uPara, 
    - tauU * (uPara * a + deltaU), 
    - tauL * (uStart * alphaL + deltaL)]

histini1= [0, 0, 0, hStart, 
    0, hPara, 
    - tauU * (hPara * a + deltaU), 
    - tauL * (hStart * alphaL + deltaL)]


u02 = [0, 0, 0, uStart, 
    0, uPara, 
    - tauL * (uPara * alphaL + deltaL + a * uPara), 
    - tauJ * (uPara * alphaL + deltaL)]

histini2 = [0, 0, 0, hStart, 
    0, hPara, 
    - tauL * (hPara * alphaL + deltaL + a * hPara), 
    - tauJ * (hPara * alphaL + deltaL)]


u03 = [0, 0, 0, uStart, 
    0, uPara, 
    - tauM * (a * uPara + deltaM), 
    - tauL * (alphaL * hStart + deltaL)]

histini3 = [0, 0, 0, hStart, 
    0, hPara, 
    - tauM * (a * hPara + deltaM), 
    - tauL * (alphaL * hStart + deltaL)]


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
para_run0 = (u00,histini0,tspan,alg, pdt)
para_run1 = (u01,histini1,tspan,alg, pdt)
para_run2 = (u02,histini2,tspan,alg, pdt)
para_run3 = (u03,histini3,tspan,alg, pdt)

# para save 
    # numéro du modèle 
    # présence de compétition chez les larves 
    # sauvgarde ou non 

para_save0 = (0,true,true) 
para_save1 = (1,true,true) 
para_save2 = (2,true,true) 
para_save3 = (3,true,true) 

# plot 
    # tuple du numéro des variables tracé
    # tuple du nom des variables 
    # taille du segment tracé
para_plot = ((2),("Larva"),500)


# Appelle de la fonction dynamique 
@time rez0 = module_analyse.traceDynamique(bc_model0, para, para_run0, para_save0, para_plot);
module_analyse.printTemp()

@time rez1 = module_analyse.traceDynamique(bc_model1, para, para_run1, para_save1, para_plot);
module_analyse.printTemp()

@time rez2 = module_analyse.traceDynamique(bc_model2, para, para_run2, para_save2, para_plot);
module_analyse.printTemp()

@time rez3 = module_analyse.traceDynamique(bc_model3, para, para_run3, para_save3, para_plot);
module_analyse.printTemp()


t_plot = 500
t_values = rez1[2]
L0 = [u[2] for u in rez0[1]]
L1 = [u[2] for u in rez1[1]]
L2 = [u[2] for u in rez2[1]]
L3 = [u[2] for u in rez3[1]]

labels = ["Without parasitoid","Eggs parasitoid","Larval parasitoid","Pupal parasitoid"]

teinte = ColorSchemes.colorschemes[:twelvebitrainbow];
positions = collect(range(0, 1; length=12))
indices = [12,6, 10, 2]
colors = [RGB{Float64}(get(teinte, positions[i])) for i in indices];
    
# Séléctionner les valeurs tracé
if t_values[end] > t_plot
    t_last_values = [x for x in t_values if x > (t_values[end] - t_plot)]
    L0 = L0[end-length(t_last_values)+1:end] 
    L1 = L1[end-length(t_last_values)+1:end] 
    L2 = L2[end-length(t_last_values)+1:end] 
    L3 = L3[end-length(t_last_values)+1:end] 
end 


t_last_values = t_last_values[1:Int(round((0.01/pdt))):end] 
L0 = L0[1:Int(round((0.01/pdt))):end]
L1 = L1[1:Int(round((0.01/pdt))):end]
L2 = L2[1:Int(round((0.01/pdt))):end]
L3 = L3[1:Int(round((0.01/pdt))):end]
u_values = [L0,L1,L2,L3]

# Créer le plot vide
plt = Plots.plot(title="Dynamics of larval density for different parasitoid", 
xlabel="Time", ylabel="Larvae density", lw=2 , grid = false);

    
# Ajouter les courbes
for (i) in 1:4
    Plots.plot!(plt,t_last_values,u_values[i],color=colors[i], linewidth = 1.5, label = labels[i])
end
display(plt)


