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

# Fonction pour résoudre le modèle et calculer les métriques
function simulate_metrics(q_values, tspan, p_fixed, u0, h, lags)
    max_vals = []
    amplitudes = []
    avg_vals = []
    cycle_durations = []
    bifurcation = []

    for q in q_values
        # Mettre à jour le paramètre "q" dans les paramètres
        p = (q, p_fixed[2:end]...)  # Remplace "q" dans la liste des paramètres

        # Résoudre le problème DDE
        prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
        sol = solve(prob, MethodOfSteps(Tsit5()))

        # On collecte les valeurs après un temps transitoire (par exemple t > 800)
        u_values = [u[1] for (t, u) in zip(sol.t, sol.u) if t > Int(tspan[2]*0.95)]
        
        # Calcul des métriques
        max_val = maximum(u_values)         # Valeur maximale
        amplitude = maximum(u_values) - minimum(u_values)  # Amplitude d'oscillation
        avg_val = mean(u_values)            # Valeur moyenne

        # Détection de la durée d'un cycle (temps entre deux pics)
        peaks = [i for i in 2:length(u_values)-1 if u_values[i-1] < u_values[i] && u_values[i] > u_values[i+1]]
        if length(peaks) > 1
            cycle_duration = sol.t[peaks[end]] - sol.t[peaks[1]]  # Durée entre le premier et le dernier pic
        else
            cycle_duration = NaN  # Si pas de cycle détecté
        end

        # Stocker les résultats pour chaque valeur de q
        push!(max_vals, max_val)
        push!(amplitudes, amplitude)
        push!(avg_vals, avg_val)
        push!(cycle_durations, cycle_duration)
    end

    return max_vals, amplitudes, avg_vals, cycle_durations
end

# Configuration du problème
tspan = (0.0, 10000)
uStart = 1
hini = -tauL * (alphaL + deltaL)
u0 = [uStart, uStart, hini]
h(p, t) = [uStart, uStart, hini]

# Liste des valeurs du paramètre "q" à analyser
q_values = 5:0.1:120  # Varie q de 5 à 100 par pas de 5

# Paramètres fixes (sans "q", qui sera remplacé dynamiquement)
p_fixed = (q, alphaL, deltaU, deltaA, deltaL, deltaM, tau1, tau2, tau3, tau4, tau5)

# Obtenir les métriques
@time max_vals, amplitudes, avg_vals, cycle_durations = simulate_metrics(q_values, tspan, p_fixed, u0, h, lags);

# Tracer les deux graphiques côte à côte
plt1 = plot(q_values, max_vals, label="Valeur maximale", xlabel="q", ylabel="Valeur max", lw=2, color=:blue)
plot!(plt1, q_values, amplitudes, label="Amplitude d'oscillation", color=:red, lw=2)
plot!(plt1, q_values, avg_vals, label="Valeur moyenne", color=:green, lw=2)
plot!(plt1, q_values, max_vals - avg_vals, label="Différence valeur moyenne - maximum", xlabel="q", ylabel="Différence", lw=2, color=:orange)

plt2 = plot(q_values, cycle_durations, label="Durée d'un cycle", xlabel="q", ylabel="Durée (temps)", lw=1, color=:green)


# Afficher les deux figures côte à côte
p = plot(plt1, plt2, layout=(1, 2), size = (1000, 400))
p
savefig(p, "/Users/Quentin/Desktop/Article1/Modele0/DD/Graphs/oscillation_q.png")  