using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using MAT
using DelimitedFiles

function bc_model(du, u, h, p, t)
    deltaU, a, q, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

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
deltaA = 1 / 12.2
deltaL = 0.0048
deltaM = 0.0150
deltaP = 1
deltaU =  0.0432
a = 1
q = 106.4/12.2
lags = [tauU, tauU + tauL + tauM, tauL + tauM]

# Fonction pour résoudre le modèle et retourner les valeurs asymptotiques
function simulate_bifurcation(deltaU_values, bc_model, tspan, p_fixed, u0, h, lags)
    results = []
    for deltaU in deltaU_values
        # Mettre à jour le paramètre "q" dans les paramètres
        p = (deltaU, p_fixed[2:end]...)  # Remplace "q" dans la liste des paramètres

        # Résoudre le problème DDE
        prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
        sol = solve(prob, MethodOfSteps(Tsit5()))

        # On collecte les valeurs asymptotiques après un temps transitoire (par ex. t > 800)
        asymptotic_values = [u[1] for (t, u) in zip(sol.t, sol.u) if t > Int(tspan[2]*0.9)]
        push!(results, asymptotic_values)
    end
    return results
end

# Configuration du problème
tspan = (0.0, 10000)
uStart = 1
hini = -tauU * (a + deltaU)
u0 = [uStart, uStart, uStart, hini]
h(p, t) = [uStart, uStart, uStart, hini]

# Liste des valeurs du paramètre "q" à analyser
deltaU_values = 0.01:0.01:1  # Varie q de 5 à 100 par pas de 5

# Paramètres fixes (sans "q", qui sera remplacé dynamiquement)
p_fixed = (deltaU, a, q, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ)

# Obtenir les résultats asymptotiques
println("Récap timer")
println("Calcules des dynamique : ")
@time results = simulate_bifurcation(deltaU_values, bc_model, tspan, p_fixed, u0, h, lags);

# Tracer le diagramme de bifurcation
plt = plot(title="Diagramme de bifurcation", 
           xlabel="moralité des oeufs (infecté ou non)", 
           ylabel="Valeur asymptotique ou amplitude")

for (i, deltaU) in enumerate(deltaU_values)
    scatter!(plt, 
             fill(deltaU, length(results[i])),  # Toutes les valeurs de q (fixe) pour les résultats correspondants
             results[i],                   # Les résultats (valeurs asymptotiques ou amplitudes)
             label=false,                  # Pas de légende
             color=:black,                  # Couleur bleue
             alpha=0.03,                    # Transparence faible (valeurs possibles : 0 à 1)
             ms=2)                         # Taille des points (marker size)
end

@time display(plt)

savefig(plt, "/Users/Quentin/Desktop/Article1/Modele1/ssDD/Graphs/bifurcation_modele1ssDD_q.png")            


