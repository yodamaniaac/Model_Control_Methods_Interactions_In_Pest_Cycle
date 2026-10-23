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
deltaA = 1 / 12.2
deltaL = 0.0048
deltaM = 0.0150
deltaP = 1
a = 1
uStart = 1

tspan = (0.0, 10000)
Ord_range = 50
Abs_range = 50

lags = [tauU, tauU + tauL + tauM, tauL + tauM]


# Définir les gamme de paramètres testé
deltaUMin = 0.001  ;  deltaUMax = 0.1
qMin = 1  ;  qMax = 70

# Créer la Matrice représentant les zone de stabilité 
deltaU_range = range(deltaUMin, deltaUMax, length=Ord_range)
q_range = range(qMin, qMax, length=Abs_range)
resultsVar = zeros(length(deltaU_range), length(q_range))
resultMean = zeros(length(deltaU_range), length(q_range))
alg = Tsit5()


for (i, q) in enumerate(q_range)
    for (j, deltaU) in enumerate(deltaU_range)
        p = (a, q, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ)
        Hini = - tauU * (a + deltaU)
        u0 = [uStart, uStart, uStart, Hini]
        h(p, t) = [uStart, uStart, uStart, Hini]        
        
        # Définit le calcule et la méthode de résolution
        prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
        alg = MethodOfSteps(Tsit5())

        # Résoudre le système d'équations différentielles à retard
        sol = solve(prob, alg)

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

p = plot(p1, p2, layout = (1, 2), size = (1000, 400))
p

savefig(p, "/Users/Quentin/Desktop/Article1/Modele1/ssDD/Graphs/limite_stab_modele1ssDD_deltaU_q.png")            

writedlm( "/Users/Quentin/Desktop/Article1/Modele1/ssDD/Tables/limite_stab_modele1ssDD_deltaU_q.csv", resultsVar, ',') 
writedlm( "/Users/Quentin/Desktop/Article1/Modele1/ssDD/Tables/densité_modele1ssDD_deltaU_q.csv", resultMean, ',') 
