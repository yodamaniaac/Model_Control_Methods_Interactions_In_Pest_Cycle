########################################
#### Test perfomance Solveur maison ####
########################################

# Charger les packages et les modules 
include("/Users/Quentin/Desktop/Article1/fonction/fc_Solveur.jl")
using .solveur_maison
using Plots, StaticArrays



# Modèle à implémenter
function modele(du, u, h, p, t)
    a, q, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ = p

    # Calcul des valeurs historiques (appel unique à chaque `h(t - ...)`)
    hist = [
        h(1,t - tauU)[4],
        h(1,t - (tauU + tauL))[4],
        h(1,t - tauL)[7],
        h(1,t - (tauU + tauL + tauM))[4],
        h(1,t - (tauL + tauM))[7],
        h(1,t - tauJ)[1],
        h(1,t - tauJ)[6],
        h(1,t - tauU)[6],
    ]

    # Modèle différentiel
    du[1] = q * u[4] - a * u[1] * u[6] - deltaU * u[1] - q * hist[1] * exp(u[7])
    du[2] = q * hist[1] * exp(u[7]) - deltaL * u[2] - q * hist[2] * exp(hist[3]) * exp(-tauL * deltaL)
    du[3] = q * hist[2] * exp(hist[3]) * exp(-tauL * deltaL) - deltaM * u[3] -
            q * hist[4] * exp(hist[5]) * exp(-tauL * deltaL) * exp(-tauM * deltaM)
    du[4] = q * hist[4] * exp(hist[5]) * exp(-tauL * deltaL) * exp(-tauM * deltaM) - deltaA * u[4]
    du[5] = a * u[1] * u[6] - deltaU * u[5] - a * hist[6] * hist[7] * exp(-deltaU * tauJ)
    du[6] = a * hist[6] * hist[7] * exp(-deltaU * tauJ) - deltaP * u[6]
    du[7] = -a * (u[6] - hist[8])
end


# Création paramètres et fixation de leur paramètres
global p, tauU, tauL, tauM, tauJ, q, a, deltaU, deltaL, deltaM, deltaA, deltaP

q = 104.6 / 12.2
deltaU = 0.8
deltaL = 0.0048
deltaM = 0.0150
deltaA = 1 / 12.2
deltaP = 1
tauU = 4.72
tauL = 16
tauM = 8.92
tauJ = tauU
a = 1
    
p = (a, q, deltaU, deltaL, deltaM, deltaA, deltaP, tauU, tauL, tauM, tauJ)

lags = [tauU, (tauU + tauL), tauL, (tauU + tauL + tauM), (tauL + tauM), tauJ, tauJ, tauU]

hini = -tauU * (a + deltaU)
uStart = 1
u0 = [uStart, uStart, uStart, uStart, uStart, uStart, hini]

pdt = 0.01
tspan = (0.0, 100.0)




# Résolution
@time t, uTot = solveur_maison.solve_dde_rk2(modele, lags, tspan, u0, pdt, p);


# Plotting
@time begin
    plt = plot(title="Dynamique du système", xlabel="Temps", ylabel="Densité des stades", lw=2)
    nom_pop_plot = ["Oeuf", "Larve", "Nymphe", "Adulte", "Oeufs infecté", "Parasitoïde"]

    for i in 1:(length(uTot[1]) - 1)
        plot!(plt, t, [u[i] for u in uTot], label=nom_pop_plot[i])
    end

    display(plt)
end


###############################
## Teste impact de la taille ##
###############################


tspanMax_values = [10, 20,50,100,200,500, 1000,2000,5000] # Les valeurs de pas de temps à tester
execution_times = Float64[]  # Vecteur pour stocker les temps d'exécution

# Boucle pour tester chaque valeur de pdt
for tspanMax in tspanMax_values
    println("Testing with tspajn = $tspanMax")
    tspan = (0.0, tspanMax)
    exec_time = @elapsed t, uTot = solveur_maison.solve_dde_rk2(modele, lags, tspan, u0, pdt)  # Mesure du temps
    push!(execution_times, exec_time)  # Stocke le temps dans le vecteur
    
end

# Affichage des résultats
println("Temps d'exécution pour chaque valeur de tspanMax :")
for (i, tspanMax) in enumerate(tspanMax_values)
    println("tspanMax = $tspanMax : time = $(execution_times[i]) seconds")
end

# Trace du graphique
plot(tspanMax_values, execution_times, lw=2, marker=:circle, xlabel="tspanMax", ylabel="Temps d'exécution (secondes)",
    title="Impact de la taille de la simulation sur le temps de simulation", legend=false, xscale=:log10, yscale=:log10)



##########################
## Teste impact des pdt ##
##########################

tspan = (0.0, 10.0)
pdt_values = [0.01, 0.001, 0.0001, 0.00001] # Les valeurs de pas de temps à tester
execution_times = Float64[]  # Vecteur pour stocker les temps d'exécution


# Boucle pour tester chaque valeur de pdt
for pdt in pdt_values
    println("Testing with pdt = $pdt")
    
    exec_time = @elapsed t, uTot = solveur_maison.solve_dde_rk2(modele, lags, tspan, u0, pdt)  # Mesure du temps
    push!(execution_times, exec_time)  # Stocke le temps dans le vecteur
    
end

# Affichage des résultats
println("Temps d'exécution pour chaque valeur de pdt :")
for (i, pdt) in enumerate(pdt_values)
    println("pdt = $pdt : time = $(execution_times[i]) seconds")
end

# Trace du graphique
plot(pdt_values, execution_times, lw=2, marker=:circle, xlabel="pdt (pas de temps)", ylabel="Temps d'exécution (secondes)",
    title="Impact du pas de temps (pdt) sur le temps de simulation", legend=false, xscale=:log10, yscale=:log10)