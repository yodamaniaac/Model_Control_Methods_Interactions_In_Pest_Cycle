using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using MAT
using DelimitedFiles

function bc_model(du, u, h, p, t)
    a, q, alphaL, deltaU, deltaL, deltaM, deltaA, deltaJ, deltaP, tauU, tauL, tauM, tauJ = p

    hist1 = h(p, t - tauU)[3]
    hist2 = h(p, t - (tauU + tauL))[3]
    hist3 = h(p, t - (tauU + tauL + tauM))[3]
    hist4 = h(p, t - tauM)[5]
    hist5 = h(p, t - tauJ)[2]
    hist6 = h(p, t - tauJ)[4]
    hist7 = h(p, t - tauL)[1]
    hist8 = h(p, t - tauM)[4]


    du[1] = q * hist1 * exp(- deltaU * tauU) - (deltaL + alphaL * u[1]) * u[1] - 
    q * hist2 * exp(- deltaU * tauU) * exp(u[5])
    du[2] = q * hist2 * exp(- deltaU * tauU) * exp(u[5]) - a * u[2] * u[4] - 
    deltaM * u[2] - q * hist3 * exp(-deltaU * tauU) * exp(hist4) * exp(u[6])
    du[3] = q * hist3 * exp(- deltaU * tauU) * exp(hist4) * exp(u[6]) - deltaA * u[3]
    du[4] = a * hist5 * hist6 * exp(- deltaJ * tauJ) - deltaP * u[4]
    du[5] = - alphaL * (u[1] - hist7)
    du[6] = - a * (u[4] - hist8) 
end
 

deltaA = 1/12.2
tauU = 4.72
tauL = 16
tauM = 8.92
tauJ = 8.92
alphaL =  0.001
deltaU =  0.0432
deltaL = 0.0048
deltaM = 0.0150
deltaP = 8
deltaJ = 0.0150
a = 1
q = 104.6/12.2

tspan = (0.0, 25.0)

uStart = 1
u0 = [uStart, uStart, uStart, uStart, 0, 0]


alg = Tsit5()
historique = 0.0
h(p, t) = [historique, historique, historique, historique, historique, historique]



lags = [tauU, tauU + tauU + tauL + tauM, tauM, tauJ, tauJ, tauL, tauM]


p = ( a, q, alphaL, deltaU, deltaL, deltaM, deltaA, deltaJ, deltaP, tauU, tauL, tauM, tauJ)


# Définit le calcule et la méthode de résolution

prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
alg = MethodOfSteps(RK4())

# Résoudre le système d'équations différentielles à retard
sol = solve(prob, alg)

# Extraire les solutions
t_values = sol.t # Les temps où la solution est évaluée
u_values = sol.u # Les valeurs des solutions correspondantes


# Noms personnalisés pour les courbes
curve_names = [ "Larves","Nymphe", "Adultes", "Parasitoïdes"]

# Créer le tracé initial
plt = plot(title="Dynamique du système (1000 derniers points)", xlabel="Temps", ylabel="Valeurs des variables", lw=2)

# Ajouter chaque composante de la solution au tracé avec les noms personnalisés
for i in 1:4
    plot!(plt, t_values, [u[i] for u in u_values], label=curve_names[i])
end

# Afficher le tracé
display(plt)
