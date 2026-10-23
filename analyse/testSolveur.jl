using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using DataFrames, CSV
using Dates
using Plots, StaticArrays, Statistics


function f(du, u, h, p, t)
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

tauU = 4.72
tauL = 16
tauM = 8.92
tau1 =  tauU
tau2 = tauU + tauL
tau3 =  tauU + tauL + tauM
retard_vect = tau3
tau4 =  tauM
tau5 =  tauL
lags = [tau1, tau2, tau3, tau4, tau5]
alphaL =  0.001
#deltaU = 0
deltaU =  0.0432
deltaL = 0.0048
deltaM = 0.0150
#deltaM = 0.53
q = 104.6/12.2
#q = 30
deltaA = 1/12.2
#deltaA = 0.75
p = (q, alphaL, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM, tau1, tau2, tau3, tau4, tau5)


hini = - tauL * (alphaL + deltaL)
uStart = 1
u0 = [0, uStart, hini]
h0= [0, 0, - tauL * deltaL]


tspan = (0.0, 60.0)
alg = "maison"
pdt = 0.01
duree_stop = 200



    t_start, t_end = tspan                                  # départ et fin de la simulation 
    n_steps = Int(ceil((t_end - t_start) / pdt) + 1)        # nb de pdt total
    retard_steps = Int(ceil(retard_vect / pdt) + 1)         # taille du retard en pas de temps 
    
    # Initialisation de U_hist et U_act
    # TODO : enlever l'ini de U_hist si on transfert juste U_act dedans 
    U_hist = Vector{SVector{length(u0), Float64}}(undef, retard_steps)      # Création d'un vecteur contenant l'historique du modèle sur la durée de la sommme des retard 
    for i in 1:retard_steps                                 # Remplir U_hist avec des SVector de zéros
        U_hist[i] = SVector{length(u0), Float64}(h0)
    end ; 
    U_hist[1]

    U_act = Vector{SVector{length(u0), Float64}}()          # Création d'un vecteur contenant les valeurs actuelle du modèle sur la durée de la sommme des retard                     
    REZ = []                                                # vecteur qui contiendra les valeur à sortir du solveur
    t = t_start                                             # valeur du temps à l'origine 
    cycle_cut = Int(ceil(duree_stop / pdt))                 # nb de pdt entre deux de vérification d'équilibre
    ratio_save = 0.01                                       # proportion de pdt à save

    eq_type = 3                                             # initialisation du type d'équilibre                                         
    
    k1 = similar(SVector{length(u0), Float64}(u0))          # initialisation des valeurs intermediaire du RK2
    k2 = similar(SVector{length(u0), Float64}(u0))
    u0

    # Fonction récupérant l'état du système pour t_retard dans le passé 
    # t_retrd : 
    # X : inutile, juste la pour homogénéiser la forme entre le solveur maison et de Julia 
    function h(X, t_retard)
        if t_retard >= t_start
            idx = Int(round((t_retard) / pdt)) + 1
            return idx <= length(U_act) ? U_act[idx] : U_hist[idx-length(U_act)]
        else
            return SVector{length(h0)}(h0)
        end
    end

    push!(U_act, SVector{length(u0)}(u0))   # integrer la valeur initial dans le vecteur U_act
    push!(REZ, u0[1])                       # integrer la valeur initial dans REZ
    
    
    for i in 1:n_steps-1                      # boucle de résolution 
        t = i * pdt                         # actualisé la valeur de t
        f(k1, U_act[end], h, p, t)          # appeller le système de DDE
        u_prelim = U_act[end] + pdt * k1    # sortir la valeur péliminaire du système 
        t_next = t + pdt                    # valeur du prochain pdt 
        f(k2, u_prelim, h, p, t_next)       # appeller le système de DDE avec initialisé pour le prochain pas de temps 
        u_new = U_act[end] + (pdt / 2) * (k1 + k2)  # trouver l'état du système après le pdt en moyennant nos deux sorties 
        push!(U_act, u_new)                 # integrer l'état du système dans U_act

        # Si U_act dépasse la taille de retard, basculer vers U_hist
        if length(U_act) > retard_steps             # si U_act plus long que la somme des retard
            U_hist = copy(U_act)                    # alors on copie indépendante U_act dans U_hist
            empty!(U_act)                           # vide U_act en gardant sa forme (mais avec une longueure de 0)
            push!(U_act, pop!(U_hist))                  # transfert la dernière valeur de U_hist à U_act
        end
    end 

    U_act 
    U_hist
    t
