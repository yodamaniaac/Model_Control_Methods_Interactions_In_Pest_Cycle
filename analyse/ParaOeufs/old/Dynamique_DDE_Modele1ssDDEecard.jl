using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using MAT
using DelimitedFiles

function bc_model(du, u, h, p, t)
    q, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

    hist1 = h(p, t - tauU)[4]
    hist2 = h(p, t - (tauU + tauL))[4]
    hist3 = h(p, t - tauL)[7] 
    hist4 = h(p, t - (tauU + tauL + tauM))[4]
    hist5 = h(p, t - (tauL + tauM))[7]
    hist6 = h(p, t - tauJ)[1]
    hist7 = h(p, t - tauJ)[6]
    hist8 = h(p, t - tauU)[6]


    du[1] = q * u[4] - a * u[1] * u[6] - deltaU * u[1] - q * hist1 * exp(u[7])
    du[2] = q * hist1 * exp(u[7]) - deltaL * u[2] - q * hist2 * exp(hist3) * exp(- tauL * deltaL)
    du[3] = q * hist2 * exp(hist3) * exp(- tauL * deltaL) - deltaM * u[3] - q * hist4 * exp(hist5) * exp(- tauL * deltaL) * exp(- tauM * deltaM)
    du[4] = q * hist4 * exp(hist5) * exp(- tauL * deltaL) * exp(- tauM * deltaM)   -   deltaA * u[4] 
    du[5] = a * u[1] * u[6] - deltaU * u[5] - a * hist6 * hist7 * exp(- deltaU * tauJ)
    du[6] = a * hist6 * hist7 * exp(- deltaU * tauJ) - deltaP * u[6]
    du[7] = - a * (u[6] - hist8) 
end
 

tauU = 4.72
tauL = 16
tauM = 8.92
tauJ = 4.72
deltaA = 1/12.2
deltaL = 0.0048
deltaM = 0.0150
deltaP = 1
a = 1
q = 106.4 / 12.2

uStart = 1
u0 = [uStart, uStart, uStart, 0]

historique = 0.0
h(p, t) = [historique, historique, historique, historique]

tspan = (0.0,10000)

lengthU = 100
MaxiMax = zeros(lengthU)
MiniMin = zeros(lengthU)
Moyenne = zeros(lengthU)

deltaU_range = range(0, 0.4, length = lengthU)

for (i,deltaU) in enumerate(deltaU_range) 

    lags = [tauU, tauU + tauL + tauM, tauL + tauM]
    p = (a, q, deltaL, deltaM, deltaA, deltaU, deltaP, tauU, tauL, tauM, tauJ)


    # Définit le calcule et la méthode de résolution

    prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
    alg = MethodOfSteps(Rodas5())

    # Résoudre le système d'équations différentielles à retard
    sol = solve(prob, alg, maxiters=1e7, abstol=1e-8, reltol=1e-6)

        # Extraire les solutions
    t_values = sol.t # Les temps où la solution est évaluée
    u_values = sol.u # Les valeurs des solutions correspondantes    

    # Prendre les 100 dernières valeurs
    last_values_indices = max(length(t_values) - 1000, 1):length(t_values)  # Assure qu'on ne dépasse pas les limites
    u_last_values = [u_values[j] for j in last_values_indices]
    u_last_values_FC = map(x -> x[1], u_last_values)
    MaxiMax[i] = maximum(u_last_values_FC)
    MiniMin[i] = minimum(u_last_values_FC)
    Moyenne[i] = mean(u_last_values_FC)
end 

# Tracé
plt = plot(deltaU_range, MaxiMax, label="MaxiMax", lw=2, color=:blue)
plot!(plt, deltaU_range, MiniMin, label="MiniMin", lw=2, color=:red)
plot!(plt, deltaU_range, Moyenne, label="Moyenne", lw=2, color=:green)

# Personnalisation
title!(plt, "Taille de la Population d'oeufs en fonction de la mortalité des oeufs")
xlabel!(plt, "Taux de Mortalité des Oeufs")
ylabel!(plt, "Taille de Pop")
plot!(plt, legend=:topright)  # Position de la légende

# Affichage
display(plt)