using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using DataFrames, CSV



function bc_model(du, u, h, p, t)
    q, alphaL, deltaU, deltaA, deltaL, deltaM, tau1, tau2, tau3, tau4, tau5 = p

    hist1 = h(p, t - tau1)[2]
    hist2 = h(p, t - tau2)[2]
    hist3 = h(p, t - tau3)[2] 
    hist4 = h(p, t - tau4)[3]
    hist5 = h(p, t - tau5)[1]

    du[1] = q * hist1 * exp(- deltaU * tauU)   -   q * hist2 * exp(-deltaU * tauU) * exp(u[3])   -   (deltaL + alphaL * u[1]) * u[1] 
    du[2] = q * hist3 * exp(- deltaU * tauU) * exp(-deltaM * tauM) * exp(hist4)   -   deltaA * u[2] 
    du[3] = - alphaL * (u[1] - hist5) 
end

tauU = 4.72
tauL = 16
tauM = 8.92
tau1 =  tauU
tau2 = tauU + tauL
tau3 =  tauU + tauL + tauM
tau4 =  tauM
tau5 =  tauL
lags = [tau1, tau2, tau3, tau4, tau5]
alphaL =  0.001
deltaU =  0.0432
deltaL = 0.0048
deltaM = 0.0150
ADD = 1

tspan = (0.0, 10000.0)
Lg_range = 3
uStart = 1
hini = - tauL * (alphaL + deltaL)
u0 = [uStart, uStart, 0]
historique = 1.0
h(p, t) = [historique, historique, 0]


# Définir les gamme de paramètres testé
deltaAMin = 0.001  ;  deltaAMax = 1
qMin = 0.001 ;  qMax = 70

# Créer la Matrice représentant les zone de stabilité 
deltaA_range = range(deltaAMin, deltaAMax, length=Lg_range)
q_range = range(qMin, qMax, length=Lg_range)
resultsVar = zeros(length(deltaA_range), length(q_range))
resultMean = zeros(length(deltaA_range), length(q_range))
alg = Tsit5()


@time for (i, q) in enumerate(q_range)
    for (j, deltaA) in enumerate(deltaA_range)
        p = (q, alphaL, deltaU, deltaA, deltaL, deltaM, tau1, tau2, tau3, tau4, tau5)

        # Définit le calcule et la méthode de résolution
        u0 = [1,1, - tauL * (alphaL + deltaL)]
        h(p,t) = u0
        prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)

        # Résolution avec un pas fixe et Runge-Kutta 2
        pdt = 0.01  # Pas de temps fixe
        sol = solve(prob,alg = MethodOfSteps(RK4()), dt=pdt)

        # Calculer la variance des variables sur les 100 derniers pas de temps
        if (length(sol[1,:]) > 100)

            resultMean[j, i] = round(mean(sol[1,end-99:end]))
            resultsVar[j, i] = std(sol[1,end-99:end]) / resultMean[j, i]

        else
            resultsVar[i, j] = 0
        end
    end
end

# Identifier les indices où resultsMean == 0
resultsVar[resultMean .== 0] .= 0


p1 = heatmap(q_range, deltaA_range, resultsVar, color=:viridis, 
xlabel="q", ylabel="deltaA", title="Carte de stabilité")
plot(p1,size = (500, 400))

p2 = heatmap(q_range, deltaA_range, resultMean, color=:viridis, 
xlabel="q", ylabel="deltaA", title="Carte des valeurs à l'Équilibre")

p = plot(p1, p2, layout = (1, 2), size = (1000, 400))
p
#savefig(p, "/Users/Quentin/Desktop/Article1/Modele0/DD/Graphs/limite_stab_modele0_deltaA_q.png")            


#writedlm( "/Users/Quentin/Desktop/Article1/Modele0/DD/Tables/limite_stab_modele0_deltaA_q.csv", resultsVar, ',') 
#writedlm( "/Users/Quentin/Desktop/Article1/Modele0/DD/Tables/densité_modele0_deltaA_q.csv", resultMean, ',') 