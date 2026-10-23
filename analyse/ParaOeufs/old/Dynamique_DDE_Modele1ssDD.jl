using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using MAT
using DelimitedFiles

function bc_model(du, u, h, p, t)
    a, q, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

    hist1 = h(p, t - tauU)[2]
    hist2 = h(p, t - (tauU + tauL + tauM))[2]
    hist3 = h(p, t - (tauL + tauM))[4]
    hist4 = h(p, t - tauJ)[1]
    hist5 = h(p, t - tauJ)[3]
    hist6 = h(p, t - tauU)[3]


    du[1] = q * u[2] - a * u[1] * u[3] - deltaU * u[1] - q * hist1 * exp(u[4])
    du[2] = q * hist2 * exp(hist3) * exp(-deltaL * tauL) * exp(-deltaM * tauM) - deltaA * u[2]
    du[3] = a * hist4 * hist5 * exp(-deltaU * tauJ) - deltaP * u[3]
    du[4] = - a * (u[3]-hist6) 
end
 

tauU = 4.72
tauL = 16
tauM = 8.92
tauJ = 4.72
deltaA = 1/12.2
deltaU = 0.0432
deltaL = 0.0048
deltaM = 0.0150
deltaP = 1
a = 1
q = 106.4 / 12.2


tspan = (0.0,10000)

hini = - tauU * (a + deltaU)
u0 = [uStart, uStart, uStart, hini]
h(p, t) = [uStart, uStart, uStart, hini]


lags = [tauU, tauU + tauL + tauM, tauL + tauM]
p = (a, q, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ)


# Définit le calcule et la méthode de résolution

prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
alg = MethodOfSteps(Tsit5())

# Résoudre le système d'équations différentielles à retard
sol = solve(prob, alg, maxiters=1e7, abstol=1e-8, reltol=1e-6)

# Extraire les solutions
t_values = sol.t # Les temps où la solution est évaluée
u_values = sol.u # Les valeurs des solutions correspondantes

# Prendre les 100 dernières valeurs
last_values_indices = max(length(t_values) - 1000, 1):length(t_values)  # Assure qu'on ne dépasse pas les limites
t_last_values = t_values[last_values_indices]
u_last_values = [u_values[i] for i in last_values_indices]

# Noms personnalisés pour les courbes
curve_names = ["Oeufs", "Adultes", "Parasitoïdes"]

# Créer le tracé initial
plt = plot(title="Dynamique du système", xlabel="Temps", ylabel="Densité des stades", lw=2)

# Ajouter chaque composante de la solution au tracé avec les noms personnalisés
for i in 1:3
    plot!(plt, t_last_values, [u[i] for u in u_last_values], label=curve_names[i])
end

# Afficher le tracé
display(plt)
