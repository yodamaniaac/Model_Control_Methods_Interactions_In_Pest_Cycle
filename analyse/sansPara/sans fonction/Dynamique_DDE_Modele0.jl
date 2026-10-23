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
 


deltaA = 1/12.2
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
q = 104.6/12.2

t_plot = 1000

tspan = (0.0, 1000)

uStart = 1
hini = - tauL * (alphaL + deltaL)
u0 = [uStart, uStart, hini]
historique = 0.0
h(p, t) = [uStart, uStart, hini]


p = (q, alphaL, deltaU, deltaA, deltaL, deltaM, tau1, tau2, tau3, tau4, tau5)


# Définit le calcule et la méthode de résolution

prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
alg = MethodOfSteps(Tsit5())

# Résoudre le système d'équations différentielles à retard
@time sol = solve(prob, alg, dtmin=1e-5, abstol=1e-8, reltol=1e-8)

# Extraire les solutions
t_values = sol.t # Les temps où la solution est évaluée
u_values = sol.u # Les valeurs des solutions correspondantes
    # Séléctionner les valeurs tracé
    last_values_indices = max(length(t_values) - t_plot, 1):length(t_values)  # Assure qu'on ne dépasse pas les limites
    t_last_values = [x for x in t_values if x > (tspan[2] - t_plot)]
    u_last_values = u_values[end-length(t_last_values)+1:end] 

# Noms personnalisés pour les courbes
curve_names = ["Larves", "Adultes"]

# Créer le tracé initial
plt = plot(title="Dynamique du système (1000 derniers points)", xlabel="Temps", ylabel="Valeurs des variables", lw=2)

# Ajouter chaque composante de la solution au tracé avec les noms personnalisés
for i in 1:1
    plot!(plt, t_values, [u[i] for u in u_values], label=curve_names[i])
end

# Afficher le tracé
display(plt)
