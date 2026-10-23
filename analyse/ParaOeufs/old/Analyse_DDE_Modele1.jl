using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using MAT
using DelimitedFiles

function bc_model(du, u, h, p, t)
    a, q, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ, tau1, tau2, tau3, tau4 = p

    hist1 = h(p, t - tau1)[3]
    hist2 = h(p, t - tau2)[3]
    hist3 = h(p, t - tau3)[3]
    hist4 = h(p, t - tauL)[6]
    hist5 = h(p, t - tau4)[6]
    hist6 = h(p, t - tauM)[5]
    hist7 = h(p, t - tauJ)[1]
    hist8 = h(p, t - tauJ)[4]
    hist9 = h(p, t - tauL)[2]
    hist10 = h(p, t - tauU)[4]


    du[1] = q * u[3] - a * u[1] * u[4] - deltaU * u[1] - q * hist1 * exp(u[6])
    du[2] = q * hist1 * exp(u[6]) - u[2] * (deltaL + alphaL * u[2]) - q * hist2 * exp(hist4) * exp(u[5])
    du[3] = q * hist3 * exp(hist5) * exp(hist6) * exp(-deltaM * tauM) - deltaA * u[3]
    du[4] = a * hist7 * hist8 * exp(-deltaU * tauJ) - deltaP * u[4]
    du[5] = - alphaL * (u[2]-hist9)
    du[6] = - a * (u[4] - hist10) 
end
 

tauU = 4.72
tauL = 16
tauM = 8.92
tauJ = 4.72
alphaL =  0.001
deltaA = 1 / 12.2
deltaL = 0.0048
deltaM = 0.0150
deltaP = 0.1
a = 1

tspan = (0.0, 10000.0)
Ord_range = 20
Abs_range = 20
uStart = 1
u0 = [uStart, uStart, uStart, uStart, 0, 0]

# Définir les gamme de paramètres testé
deltaUMin = 0.001  ;  deltaUMax = 0.1
qMin = 1  ;  qMax = 70

# Créer la Matrice représentant les zone de stabilité 
deltaU_range = range(deltaUMin, deltaUMax, length=Ord_range)
q_range = range(qMin, qMax, length=Abs_range)
resultsVar = zeros(length(deltaU_range), length(q_range))
resultMean = zeros(length(deltaU_range), length(q_range))
alg = Tsit5()

historique = 0.0
h(p, t) = [historique, historique, historique, historique, historique, historique]


for (i, q) in enumerate(q_range)
    for (j, deltaU) in enumerate(deltaU_range)
        tau1 = tauU
        tau2 = tauU + tauL
        tau3 = tauU + tauL + tauM
        tau4 = tauL + tauM
        lags = [tau1, tau2, tau3, tau4]
        p = (a, q, alphaL, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ, tau1, tau2, tau3, tau4)


        # Définit le calcule et la méthode de résolution

        prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
        alg = MethodOfSteps(Tsit5())

        # Résoudre le système d'équations différentielles à retard
        sol = solve(prob, alg,reltol=1e-6, abstol=1e-6)

        # Calculer la variance des variables sur les 100 derniers pas de temps
        if (length(sol[1,:]) > 100)
            
            resultMean[j, i] = round(mean(sol[1,end-99:end]))
            resultsVar[j, i] = std(sol[1,end-99:end]) / resultMean[j, i]
            
            if resultsVar[j, i] < 0
                resultsVar[j, i] = 0
            end

        else
            resultsVar[i, j] = 0
        end
    end
end

# Pour sauvgarder la matrice
# writedlm("/Users/Quentin/Desktop/codes/test julia/Analyse_DDE/matrix_Var.csv", resultsVar, ',') 
# writedlm("/Users/Quentin/Desktop/codes/test julia/Analyse_DDE/matrix_Mean.csv", resultMean, ',') 
p1 = heatmap(q_range, deltaU_range,resultsVar, color=:viridis, 
xlabel="q", ylabel="deltaU", title="Carte de stabilité")

p2 = heatmap(q_range, deltaU_range,resultMean, color=:viridis, 
xlabel="q", ylabel="deltaU", title="Carte des valeurs à l'Équilibre")

plot(p1, p2, layout = (1, 2), size = (1000, 400))
