#####################################################
#### Analyse Modèle Ravageur 4 stades (Modèle 0) ####
#####################################################


# Import packages
if Sys.KERNEL == :NT
    include("P:/Article1/fonction/fc_Solveur.jl")
    include("P:/Article1/fonction/fc_dynamique.jl")
else
    include("/Users/Quentin/Desktop/Article1/fonction/fc_Solveur.jl")
    include("/Users/Quentin/Desktop/Article1/fonction/fc_dynamique.jl")
end

using .solveur_maison
using .module_analyse
using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using DataFrames, CSV
using Dates
using FFTW


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
q = 30
deltaA = 0.25


p = (q, alphaL, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM, tau1, tau2, tau3, tau4, tau5)


hini = - tauL * (alphaL + deltaL)
uStart = 1
u0 = [uStart, uStart, hini]
histini= [0, 0, - tauL * alphaL]

tspan = (0.0, 20000.0)
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
para_run = (u0,histini,tspan,"maison",pdt)

# para save 
    # numéro du modèle 
    # présence de compétition chez les larves 
    # sauvgarde ou non 
para_save = (0,true,false) 

# plot 
    # tuple du numéro des variables tracé
    # tuple du nom des variables 
    # taille du segment tracé
para_plot = ((1),("Larve","Adulte"),10000)


# Appelle de la fonction dynamique 
@time teste = module_analyse.traceDynamique(bc_model, para, para_run, para_save, para_plot)

# Exemple : données temporelles pour une variable
t = teste[2][Int(round(end/2+1)):Int(round((1/pdt))):end]                   # Temps (issu de votre simulation)
u = [u[1] for u in teste[1]][Int(round(end/2+1)):Int(round((1/pdt))):end]      # Valeurs pour la première variable
plot(t[end-1000:end],u[end-1000:end])

# Vérifier l'uniformité des pas de temps 
dt = t[2] - t[1]             # Pas de temps (fixe)
is_uniform = all(x -> x ≈ dt, diff(t))  


# Transformée de Fourier
U_fft = fft(u)  # Transformée
frequencies = fftfreq(length(U_fft), dt)  # Fréquences associées
amplitudes = abs.(U_fft)  # Amplitudes (module des coefficients)

# Trouver la fréquence max 
dominant_frequency = frequencies[argmax(amplitudes[2:end÷2])]
dominant_period = 1/dominant_frequency

# Visualisation du spectre
plot(frequencies[1:end÷2], amplitudes[1:end÷2],  # Partie positive du spectre
     xlabel="Fréquence (Hz)", ylabel="Amplitude", title="Analyse spectrale")

