using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using MAT
using DelimitedFiles

function bc_model(du, u, h, p, t)
    q, alphaL, deltaU, deltaA, deltaL, deltaM, tau1, tau2, tau3, tau4, tau5 = p

    hist1 = h(p, t - tau1)[2]
    hist2 = h(p, t - tau2)[2]
    hist3 = h(p, t - tau3)[2] 
    hist4 = h(p, t - tau4)[3]
    hist5 = h(p, t - tau5)[1]

    du[1] = q  * hist1 * exp(-deltaU * tau1)   -   q  * hist2 * exp(-deltaU * tau1) * exp(u[3])   -   (deltaL + alphaL * u[1]) * u[1] 
    du[2] = q  * hist3 * exp(-deltaU * tau1) * exp(-deltaM * tau4) * exp(hist4)   -   deltaA * u[2] 
    du[3] = alphaL * (-u[1] + hist5) 
end
 


deltaA = 1
tauU = 4.72
tauL = 16
tauN = 8.92
tau1 =  tauU
tau2 = tauU + tauL
tau3 =  tauU + tauL + tauN
tau4 =  tauN
tau5 =  tauL
lags = [tau1, tau2, tau3, tau4, tau5]
alphaL =  0.001
deltaU =  0.0432
deltaL = 0.0048
deltaM = 0.0150
q= 106.4

# Fonction pour résoudre le modèle et retourner les valeurs asymptotiques
function simulate_bifurcation(q_values, tspan, p_fixed, u0, h, lags)
    results = []
    for q in q_values
        # Mettre à jour le paramètre "q" dans les paramètres
        p = (q, p_fixed[2:end]...)  # Remplace "q" dans la liste des paramètres

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
hini = -tauL * (alphaL + deltaL)
u0 = [uStart, uStart, hini]
h(p, t) = [uStart, uStart, hini]

# Liste des valeurs du paramètre "q" à analyser
q_values = 5:0.1:100  # Varie q de 5 à 100 par pas de 5

# Paramètres fixes (sans "q", qui sera remplacé dynamiquement)
p_fixed = (q,alphaL, deltaU, deltaA, deltaL, deltaM, tau1, tau2, tau3, tau4, tau5)

# Obtenir les résultats asymptotiques
println("Récap timer")
println("Calcules des dynamique : ")
@time results = simulate_bifurcation(q_values, tspan, p_fixed, u0, h, lags);

# Tracer le diagramme de bifurcation
plt = plot(title="Diagramme de bifurcation", 
           xlabel="q (Taux de production)", 
           ylabel="Valeur asymptotique ou amplitude")

for (i, q) in enumerate(q_values)
    scatter!(plt, 
             fill(q, length(results[i])),  # Toutes les valeurs de q (fixe) pour les résultats correspondants
             results[i],                   # Les résultats (valeurs asymptotiques ou amplitudes)
             label=false,                  # Pas de légende
             color=:black,                  # Couleur bleue
             alpha=0.03,                    # Transparence faible (valeurs possibles : 0 à 1)
             ms=2)                         # Taille des points (marker size)
end

@time display(plt)

savefig(plt, "/Users/Quentin/Desktop/Article1/Modele0/DD/Graphs/bifurcation_modele0_q.png")            


# enquête sur les zones grise et les zones de traits noir 
p1 = (79, alphaL, deltaU, deltaA, deltaL, deltaM, tau1, tau2, tau3, tau4, tau5)
p2 = (81, alphaL, deltaU, deltaA, deltaL, deltaM, tau1, tau2, tau3, tau4, tau5)

prob1 = DDEProblem(bc_model, u0, h, tspan, p1; constant_lags = lags)
prob2 = DDEProblem(bc_model, u0, h, tspan, p2; constant_lags = lags)
alg = MethodOfSteps(Tsit5())

# Résoudre le système d'équations différentielles à retard
@time sol1 = solve(prob1, alg, dtmin=1e-5, abstol=1e-8, reltol=1e-8);
@time sol2 = solve(prob2, alg, dtmin=1e-5, abstol=1e-8, reltol=1e-8);

# Extraire les solutions
t_values1 = sol1.t; # Les temps où la solution est évaluée
u_values1 = sol1.u; # Les valeurs des solutions correspondantes
t_values1 = t_values1[end-1099:end];
u_values1 = u_values1[end-1099:end];

t_values2 = sol2.t; # Les temps où la solution est évaluée
u_values2 = sol2.u; # Les valeurs des solutions correspondantes
t_values2 = t_values2[end-1099:end];
u_values2 = u_values2[end-1099:end];

# Noms personnalisés pour les courbes
curve_names1 = ["Larves Gris", "Adultes Gris"]
curve_names2 = ["Larves Trait", "Adultes Trait"]

# Créer le tracé initial
plt = plot(title="Comparaison dynamique zonne grise / zone noire", xlabel="Temps", ylabel="Valeurs des variables", lw=2)

# Ajouter chaque composante de la solution au tracé avec les noms personnalisés
for i in 1:2
    p1=plot!(plt, t_values1, [u[i] for u in u_values1], label=curve_names1[i])
    p2=plot!(plt, t_values2, [u[i] for u in u_values2], label=curve_names2[i])
end

# Afficher le tracé
display(plt)
