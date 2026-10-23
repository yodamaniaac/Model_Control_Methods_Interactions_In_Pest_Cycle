#######################################
#### Répertoir des Solveurs Maison ####
#######################################

### Code générale ---- 
module solveur_maison

# Fonction exporté 
export solve_dde_rk2, solve_dde_rk2_explicit_mid_point, solve_dde_rk2_first,solve_dde_rk2_Chat

# Import
using Plots, StaticArrays, Statistics


### Création de la fonction Solveur RK2 ----

# Paramètres 
#   f   :  fonction du modèle
#
#   retard_vect : vecteur des retards
#
#   tspan : tuple de 2 Float64: 
#       tps au départ de simul
#       tps de fin de simul
#
#   u0  :  ajout initiale dans les pops 
#
#   pdt : Float64 du pas de temps utilisé
function solve_dde_rk2_explicit_mid_point(f, retard_vect, tspan, u0, pdt, p)
    
    t_start, t_end = tspan
    n_steps = Int(round((t_end - t_start) / pdt) + 1)
    t = range(t_start, stop=t_end, length=n_steps)
    uTot = Vector{SVector{length(u0), Float64}}(undef, n_steps)
    uTot[1] = SVector{length(u0)}(u0)
    # Préallocation des variables temporaires
    k1 = similar(uTot[1])
    k2 = similar(uTot[1])

    # Fonction d'accès aux valeurs retardées
    function h(X, t_retard)
        if t_retard >= t_start
            idx = Int(round((t_retard - t_start) / pdt)) + 1
            return uTot[clamp(idx, 1, n_steps)]
        else
            return SVector{length(u0)}(u0)
        end
    end

    for i in 1:n_steps - 1

        f(k1, uTot[i], h, p, t[i])
        uTot_mid = uTot[i] + 0.5 * pdt * k1
        t_mid = t[i] + 0.5 * pdt
        f(k2, uTot_mid, h, p, t_mid)
        uTot[i + 1] = uTot[i] + pdt * k2
    end

    return t, uTot
end




### Création de la fonction Solveur RK2 ----

# Paramètres 
#   f   :  fonction du modèle
#
#   retard_vect : float de la somme des vecteur des retards
#
#   tspan : tuple de 2 Float64: 
#       tps au départ de simul
#       tps de fin de simul
#
#   u0  :  ajout initiale dans les pops
#
#   pdt : Float64 du pas de temps utilisé
function solve_dde_rk2_first(f, retard_vect, tspan, u0, pdt, p, h0, duree_stop = 200)
    t_start, t_end = tspan
    n_steps = Int(round((t_end - t_start) / pdt) + 1)
    t = range(t_start, stop=t_end, length=n_steps)
    uTot = Vector{SVector{length(u0), Float64}}(undef, n_steps)
    uTot[1] = SVector{length(u0)}(u0)
    eq_type = 3

    # Préallocation des variables temporaires
    k1 = similar(uTot[1])
    k2 = similar(uTot[1])
    
    # Fonction d'accès aux valeurs retardées
    function h(X, t_retard)
        if t_retard >= t_start
            idx = Int(round((t_retard - t_start) / pdt)) + 1
            return uTot[clamp(idx, 1, n_steps)]
        else
            return SVector{length(h0)}(h0)
        end
    end

    for i in 1:n_steps - 1
        # Resolution RK2 Heun
        f(k1, uTot[i], h, p, t[i])
        u_prelim = uTot[i] + pdt * k1
        t_next = t[i] + pdt
        
        f(k2, u_prelim, h, p, t_next)
        uTot[i + 1] = uTot[i] + (pdt / 2) * (k1 + k2) 


        # Condition d'arret
        cycle_cut = Int(round(duree_stop / pdt)) 

        if i > cycle_cut * 2.5 && i % cycle_cut == 0 

            # Condition d'arret pour equilibre stable 
            dyna_stop = [v[1] for v in uTot[i-cycle_cut:i]]

            if std(dyna_stop[Int(round(length(dyna_stop)/2)):end]) < 10^(-4)

                uTot = resize!(uTot, i + 1)
                t = t[1:i+1]
                
                if mean(dyna_stop)<10^-2
                    eq_type = 0
                    sortir_cdt = (mean(dyna_stop), std(dyna_stop), t[end],"0")
                else 
                    eq_type = 1
                    sortir_cdt = (mean(dyna_stop),std(dyna_stop), t[end],"1")
                end 
                
                println(sortir_cdt)

                break

            # Condition d'arret pour equilibre instable 

            elseif  abs(maximum(dyna_stop[Int(round(length(dyna_stop)/2)):end]) - maximum(dyna_stop[1:Int(round(length(dyna_stop)/2))-1])) < 10^-4 &&
                abs(minimum(dyna_stop[Int(round(length(dyna_stop)/2)):end]) - minimum(dyna_stop[1:Int(round(length(dyna_stop)/2))-1])) < 10^-4  

                uTot = resize!(uTot, i + 1)
                t = t[1:i+1]
                sortir_cdt = (mean(dyna_stop),std(dyna_stop), t[end],"2")
                eq_type = 2
                println(sortir_cdt)
           
                break
            end
        end 
    end
    if eq_type == 3
        cycle_cut = Int(round(duree_stop / pdt)) 
        if cycle_cut <= length(uTot)
            dyna_stop = [v[1] for v in uTot[end-cycle_cut:end]]
            sortir_cdt = (mean(dyna_stop),std(dyna_stop), t[end],"3")
        else 
            dyna_stop = [v[1] for v in uTot[1:end]]
            sortir_cdt = (mean(dyna_stop),std(dyna_stop), t[end],"3")
        end

    end 
    return t, uTot, eq_type
end








### Création de la fonction Solveur RK2 ----

# Paramètres 
#   f   :  fonction du modèle
#
#   retard_vect : float de la somme des vecteur des retards
#
#   tspan : tuple de 2 Float64: 
#       tps au départ de simul
#       tps de fin de simul
#
#   u0  :  ajout initiale dans les pops 
#
#   pdt : Float64 du pas de temps utilisé

function solve_dde_rk2(f, retard_vect, tspan, u0, pdt, p, h0, duree_stop = 200)
    t_start, t_end = tspan                                  # départ et fin de la simulation 
    n_steps = Int(ceil((t_end - t_start) / pdt) + 1)        # nb de pdt total
    retard_steps = Int(ceil(retard_vect / pdt) + 1)         # taille du retard en pas de temps 
    
    # Initialisation de U_hist et U_act
    # TODO : enlever l'ini de U_hist si on transfert juste U_act dedans 
    U_hist = Vector{SVector{length(u0), Float64}}(undef, retard_steps)      # Création d'un vecteur contenant l'historique du modèle sur la durée de la sommme des retard 
    for i in 1:retard_steps                                 # Remplir U_hist avec des SVector de zéros
        U_hist[i] = SVector{length(u0), Float64}(h0)
    end ; 

    
    U_act = Vector{SVector{length(u0), Float64}}()          # Création d'un vecteur contenant les valeurs actuelle du modèle sur la durée de la sommme des retard                     
    REZ = []                                                # vecteur qui contiendra les valeur à sortir du solveur
    t = t_start                                             # valeur du temps à l'origine 
    cycle_cut = Int(ceil(duree_stop / pdt))                 # nb de pdt entre deux de vérification d'équilibre
    ratio_save = 0.01                                       # proportion de pdt à save

    eq_type = 3                                             # initialisation du type d'équilibre                                         
    
    k1 = similar(SVector{length(u0), Float64}(u0))          # initialisation des valeurs intermediaire du RK2
    k2 = similar(SVector{length(u0), Float64}(u0))


    # Fonction récupérant l'état du système pour t_retard dans le passé 
    # t_retrd : 
    # X : inutile, juste la pour homogénéiser la forme entre le solveur maison et de Julia 
    function h(X, t_retard)
        if t_retard >= t_start                              # si le retard arrive après le départ
            idx = Int(round((t-t_retard) / pdt))            # taille du retard (non relatif à t)
            return idx < length(U_act) ? U_act[end-idx] : U_hist[end-(idx-length(U_act))]
        else
            return SVector{length(h0)}(h0)
        end
    end

    push!(U_act, SVector{length(u0)}(u0))   # integrer la valeur initial dans le vecteur U_act
    push!(REZ, u0)                       # integrer la valeur initial dans REZ
    
    
    for i in 1:n_steps-1                      # boucle de résolution 
        t = i * pdt                         # actualisé la valeur de t
        f(k1, U_act[end], h, p, t)          # appeller le système de DDE
        u_prelim = U_act[end] + pdt * k1    # sortir la valeur péliminaire du système 
        t_next = t + pdt                    # valeur du prochain pdt 
        f(k2, u_prelim, h, p, t_next)       # appeller le système de DDE avec initialisé pour le prochain pas de temps 
        u_new = U_act[end] + (pdt / 2) * (k1 + k2)  # trouver l'état du système après le pdt en moyennant nos deux sorties 
        push!(U_act, u_new)                 # integrer l'état du système dans U_act

        # Si U_act dépasse la taille de retard, basculer vers U_hist
        if length(U_act) > retard_steps             # si U_act plus long que la somme des retard
            U_hist = copy(U_act)                    # alors on copie indépendante U_act dans U_hist
            empty!(U_act)                           # vide U_act en gardant sa forme (mais avec une longueure de 0)
            push!(U_act, pop!(U_hist))                  # transfert la dernière valeur de U_hist à U_act
        end

        # Sauvegarde pour analyse
        if i % (ratio_save/pdt) == 0    
            push!(REZ, u_new)
        end

        # Condition d'arrêt
        if i > cycle_cut * 2.5 && i % cycle_cut == 0
            dyna_stop = [v[1] for v in REZ[Int(round(end - duree_stop / ratio_save)):end]]
            if std(dyna_stop[Int(round(length(dyna_stop) / 2)):end]) < 10^(-3)
                eq_type = if mean(dyna_stop) < 10^-2 0 else 1 end
                println((mean(dyna_stop), std(dyna_stop), t, string(eq_type)))
                break
            elseif abs(maximum(dyna_stop[Int(round(length(dyna_stop) / 2)):end]) - maximum(dyna_stop[1:Int(round(length(dyna_stop) / 2))-1])) < 0.0001 &&
                    abs(minimum(dyna_stop[Int(round(length(dyna_stop) / 2)):end]) - minimum(dyna_stop[1:Int(round(length(dyna_stop) / 2))-1])) < 0.0001
                eq_type = 2
                println((mean(dyna_stop), std(dyna_stop), t, "2"))
                break
            end
        end
    end

    if eq_type == 3
        dyna_stop = [v[1] for v in REZ[max(1, end - cycle_cut):end]]
        println((mean(dyna_stop), std(dyna_stop), t, "3"))
    end


    return 0:ratio_save:t, REZ, eq_type
end








### Création de la fonction Solveur RK2 ----

# Paramètres 
#   f   :  fonction du modèle
#
#   retard_vect : float de la somme des vecteur des retards
#
#   tspan : tuple de 2 Float64: 
#       tps au départ de simul
#       tps de fin de simul
#
#   u0  :  ajout initiale dans les pops 
#
#   pdt : Float64 du pas de temps utilisé
function solve_dde_rk2_second(f, retard_vect, tspan, u0, pdt, p, h0, duree_stop = 200)
    t_start, t_end = tspan
    n_steps = Int(round((t_end - t_start) / pdt) + 1)
    
    # Initialisation de U_hist et U_act
    U_hist = Vector{SVector{length(u0), Float64}}(undef, Int(round(retard_vect / pdt)))
    U_act = Vector{SVector{length(u0), Float64}}(undef, Int(round(retard_vect / pdt)))
    REZ = []
    t = 0
    # Initialisation de uTot, de l'état initial
    uTot = Vector{SVector{length(u0), Float64}}(undef, n_steps)
    uTot[1] = SVector{length(u0)}(u0)

    eq_type = 3
    k1 = similar(uTot[1])
    k2 = similar(uTot[1])

    function h(X, t_retard)
        if t_retard >= t_start
            idx = Int(round((t_retard - t_start) / pdt)) + 1
            return uTot[clamp(idx, 1, n_steps)]
        else
            return SVector{length(h0)}(h0)
        end
    end

    for i in 1:n_steps - 1
        # Résolution RK2 Heun
        t = i * pdt
        f(k1, uTot[i], h, p, t)
        u_prelim = uTot[i] + pdt * k1
        t_next = t + pdt
        f(k2, u_prelim, h, p, t_next)
        uTot[i + 1] = uTot[i] + (pdt / 2) * (k1 + k2)

        # Ajouter la nouvelle valeur à U_act
        push!(U_act, uTot[i + 1])
        
        # Si U_act est plein, remplacer U_hist et réinitialiser U_act
        if length(U_act) >= length(U_hist)
            U_hist = copy(U_act)
            empty!(U_act)
        end
        
        # Sauvegarder l'état actuel dans REZ chaque 1/pdt
        if i % Int(round(1 / pdt)) == 0
            push!(REZ, uTot[i + 1])
        end
        
        # Condition d'arrêt
        cycle_cut = Int(round(duree_stop / pdt))
        if i > cycle_cut * 2.5 && i % cycle_cut == 0
            dyna_stop = [v[1] for v in uTot[i - cycle_cut:i]]
            if std(dyna_stop[Int(round(length(dyna_stop) / 2)):end]) < 10^(-4)
                uTot = resize!(uTot, i + 1)
                t_out = 1:t
                if mean(dyna_stop) < 10^-2
                    eq_type = 0
                    sortir_cdt = (mean(dyna_stop), std(dyna_stop), t, "0")
                else
                    eq_type = 1
                    sortir_cdt = (mean(dyna_stop), std(dyna_stop), t, "1")
                end
                println(sortir_cdt)
                break
            elseif abs(maximum(dyna_stop[Int(round(length(dyna_stop) / 2)):end]) - maximum(dyna_stop[1:Int(round(length(dyna_stop) / 2))])) < 10^-4 &&
                    abs(minimum(dyna_stop[Int(round(length(dyna_stop) / 2)):end]) - minimum(dyna_stop[1:Int(round(length(dyna_stop) / 2))])) < 10^-4
                uTot = resize!(uTot, i + 1)
                t_out = 1:t
                sortir_cdt = (mean(dyna_stop), std(dyna_stop), t, "2")
                eq_type = 2
                println(sortir_cdt)
                break
            end
        end
    end

    if eq_type == 3
        cycle_cut = Int(round(duree_stop / pdt))
        if cycle_cut <= length(uTot)
            dyna_stop = [v[1] for v in uTot[end - cycle_cut:end]]
            sortir_cdt = (mean(dyna_stop), std(dyna_stop), t, "3")
        else
            dyna_stop = [v[1] for v in uTot[1:end]]
            sortir_cdt = (mean(dyna_stop), std(dyna_stop), t, "3")
        end
    end

    return 1:t, REZ, eq_type
end









# Solveur Chat
function solve_dde_rk2_Chat(f, retard_vect, tspan, u0, pdt, p, h0, duree_stop=200)
    t_start, t_end = tspan
    n_steps = Int(round((t_end - t_start) / pdt) + 1)
    t = range(t_start, stop=t_end, length=n_steps)
    uTot = Vector{SVector{length(u0), Float64}}(undef, n_steps)
    uTot[1] = SVector{length(u0)}(u0)
    eq_type = 3

    # Préallocation des variables temporaires
    k1 = similar(uTot[1])
    k2 = similar(uTot[1])
    
    # Fonction d'accès aux valeurs retardées avec interpolation linéaire
    function h(X, t_retard)
        if t_retard >= t_start
            idx = (t_retard - t_start) / pdt
            i_low = floor(Int, idx) + 1
            i_high = clamp(i_low + 1, 1, n_steps)
            alpha = idx - (i_low - 1)
            return (1 - alpha) * uTot[clamp(i_low, 1, n_steps)] + alpha * uTot[i_high]
        else
            return SVector{length(h0)}(h0)
        end
    end

    for i in 1:n_steps - 1
        # Méthode de Heun (RK2 amélioré)
        f(k1, uTot[i], h, p, t[i])
        u_prelim = uTot[i] + pdt * k1
        t_next = t[i] + pdt
        
        f(k2, u_prelim, h, p, t_next)
        uTot[i + 1] = uTot[i] + (pdt / 2) * (k1 + k2) 

        # Condition d'arrêt
        cycle_cut = Int(round(duree_stop / pdt)) 

        if i > cycle_cut * 2.5 && i % cycle_cut == 0 
            dyna_stop = [v[1] for v in uTot[i-cycle_cut:i]]
            
            # Vérification stabilité (convergence absolue)
            if maximum(abs, diff(dyna_stop[Int(round(length(dyna_stop)/2)):end])) < 10^(-4)
                uTot = resize!(uTot, i + 1)
                t = t[1:i+1]
                if mean(dyna_stop) < 10^-2
                    eq_type = 0 
                else 
                    eq_type = 1
                end
                println((mean(dyna_stop), std(dyna_stop), t[end], string(eq_type)))
                break
            end
            
            # Vérification cycle limite avec autocorrélation
            autocorr = sum((dyna_stop .- mean(dyna_stop)) .* (circshift(dyna_stop, cycle_cut ÷ 2) .- mean(dyna_stop)))
            if abs(autocorr) > 0.9 * sum(abs2, dyna_stop)
                uTot = resize!(uTot, i + 1)
                t = t[1:i+1]
                println((mean(dyna_stop), std(dyna_stop), t[end], "2"))
                eq_type = 2
                break
            end
        end
    end

    # Si aucun critère n'a stoppé la simulation
    if eq_type == 3
        cycle_cut = Int(round(duree_stop / pdt))
        dyna_stop = [v[1] for v in uTot[max(1, end-cycle_cut):end]]
        println((mean(dyna_stop), std(dyna_stop), t[end], "3"))
    end

    return t, uTot, eq_type
end



end


