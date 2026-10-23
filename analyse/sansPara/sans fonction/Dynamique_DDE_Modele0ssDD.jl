using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using MAT
using DelimitedFiles

function bc_model(du, u, h, p, t)
    q, a, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM= p

    hist1 = h(p, t - (tauU + tauL + tauM))[1]

    du[1] = q * hist1 * exp(-deltaU * tauU) * exp(-deltaL * tauL) * exp(-deltaM * tauM) - deltaA * u[1]
end
 

tauU = 4.72
tauL = 16
tauM = 8.92
deltaA = 1/12.2
deltaU =  0.0432
deltaL = 0.0048
deltaM = 0.0150

a = 1
q = 106.4

tspan = (0.0, 1000.0)

uStart = 1
u0 = [uStart]


alg = Tsit5()
historique = 0.0
h(p, t) = [historique]



p = (q, a, deltaU, deltaL, deltaM, deltaA, tauU, tauL, tauM)
lags = [tauU + tauL + tauM]

# Définit le calcule et la méthode de résolution

prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
alg = MethodOfSteps(RK4())

# Résoudre le système d'équations différentielles à retard
sol = solve(prob, alg)

# Extraire les solutions
t_values = sol.t # Les temps où la solution est évaluée
u_values = sol.u # Les valeurs des solutions correspondantes


# Noms personnalisés pour les courbes
curve_names = ["Oeufs", "Larves", "Adultes", "Parasitoïdes"]

# Créer le tracé initial
plt = plot(title="Dynamique du système (1000 derniers points)", xlabel="Temps", ylabel="Valeurs des variables", lw=2)

# Ajouter chaque composante de la solution au tracé avec les noms personnalisés
for i in 1:1
    plot!(plt, t_values, [u[i] for u in u_values], label=curve_names[i])
end

# Afficher le tracé
display(plt)
