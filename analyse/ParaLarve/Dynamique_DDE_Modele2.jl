using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using MAT
using DelimitedFiles

function bc_model(du, u, h, p, t)
    q, a, alphaL, alphaJ, deltaU, deltaL, deltaM, deltaA, deltaJ, deltaP, tauU, tauL, tauM, tauJ = p

    hist1 = h(p, t - tauU)[2]
    hist2 = h(p, t - (tauU + tauL))[2]
    hist3 = h(p, t - (tauU + tauL + tauM))[2]
    hist4 = h(p, t - tauM)[6]
    hist5 = h(p, t - tauJ)[1]
    hist6 = h(p, t - tauJ)[3]
    hist7 = h(p, t - tauL)[4]
    hist8 = h(p, t - tauL)[1]
    hist9 = h(p, t - tauL)[3]


    du[1] = q * hist1 - a * u[1] * u[4] - (deltaL + alphaL * (u[1]+u[3])) * u[1] - q * hist2 * exp(- deltaU * tauU) * exp(u[5])
    du[2] = q * hist3 * exp(- deltaU * tauU) * exp(hist4)* exp(- deltaM * tauM) - deltaA * u[2]
    du[3] = a * u[1] * u[3] - (deltaJ + alphaJ * (u[1]+u[3])) * u[3] - a * hist5 * hist6 * exp(u[6])
    du[4] = a * hist5 * hist6 * exp(u[6]) - deltaP * u[4]
    du[5] = - a * (u[4] - hist7) - alphaL * ((u[1] + u[3]) - (hist8 + hist9))
    du[6] = - alphaJ * ((u[1] + u[3]) - (hist8 + hist9)) 
end
 

tauU = 4.72
tauL = 16
tauM = 8.92
tauJ = 16
deltaA = 1/12.2
alphaL =  0.001
alphaJ =  0.001
deltaU =  0.0432
deltaL = 0.0048
deltaM = 0.0150
deltaP = 8
deltaJ = 0.0048
a = 1
q = 106.4 /12.2

tspan = (0.0, 1000.0)

uStart = 1
u0 = [uStart, uStart, uStart, uStart, 0, 0]


alg = Tsit5()
historique = 0.0
h(p, t) = [historique, historique, historique, historique, historique, historique]


lags = [tauU, tauU + tauL, tauU + tauL + tauM, tauM, tauJ, tauJ, tauL, tauL, tauL]
p = (q, a, alphaL, alphaJ, deltaU, deltaL, deltaM, deltaA, deltaJ, deltaP, tauU, tauL, tauM, tauJ)


# Définit le calcule et la méthode de résolution

prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
alg = MethodOfSteps(RK4())

# Résoudre le système d'équations différentielles à retard
sol = solve(prob, alg)

# Extraire les solutions
t_values = sol.t # Les temps où la solution est évaluée
u_values = sol.u # Les valeurs des solutions correspondantes

# Prendre les 100 dernières valeurs
last_values_indices = max(length(t_values) - 1000, 1):length(t_values)  # Assure qu'on ne dépasse pas les limites
t_last_values = t_values[last_values_indices]
u_last_values = [u_values[i] for i in last_values_indices]

# Noms personnalisés pour les courbes
curve_names = [ "Larves", "Adultes","Larve_Parasite", "Parasitoïdes"]

# Créer le tracé initial
plt = plot(title="Dynamique du système (1000 derniers points)", xlabel="Temps", ylabel="Valeurs des variables", lw=2)

# Ajouter chaque composante de la solution au tracé avec les noms personnalisés
for i in 1:4
    plot!(plt, t_last_values, [u[i] for u in u_last_values], label=curve_names[i])
end

# Afficher le tracé
display(plt)
u_last_values