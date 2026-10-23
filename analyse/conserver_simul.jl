    # Résoudre le problème selon le solveur choisis 
    if alg == "maison"
        # Définir le problème et le résoudre
        t_values, u_values, eq_type = solveur_maison.solve_dde_rk2_first(modele, totLag, tspan, u0, pdt, p, h0)
    
    elseif alg == "fixe"
        # Définir le problème 
        frq_save = 0.001
        prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)

        # Définir une condition d'arrêt basée sur l'écart type
        function condition(u, t, integrator)
            # Vérifier toutes les 200 unités de temps
            if t >= 500.0*20000000000.5 && t % 500.0 < integrator.dt
                # Calculer l'écart type des valeurs de u durant les 200 dernières unités de temps

                # Filtrez les valeurs qui sont supérieures ou égales au seuil
                recent_time = integrator.sol.t[integrator.sol.t .>= integrator.sol.t[end]-200]
                recent_values = [v[1] for v in integrator.sol.u[end-(length(recent_time)-1):end]]
                #println("time : ", recent_time[1],"  -  ", recent_time[end])
                #println("valeur 1/end : ", recent_values[1],"  -  ", recent_values[end])
                #println("time length : ",length(recent_time))
                #println("valeur length : ",length(recent_values))
                std_u = std(recent_values)
                #println("max-max : ",(abs(maximum(recent_values[Int(round(length(recent_values)/2)):end]) - 
                #maximum(recent_values[1:Int(round(length(recent_values)/2))-1]))) )
                #println("moyenne : ", mean(recent_values) )
                if std_u < 0.001
                    println(1)
                    return true
#                elseif (abs(maximum(recent_values[Int(round(length(recent_values)/2)):end]) - 
#                    maximum(recent_values[1:Int(round(length(recent_values)/2))-1]))) / mean(recent_values) < 0.001 &&
#                    abs(minimum(recent_values[Int(round(length(recent_values)/2)):end]) - 
#                    minimum(recent_values[1:Int(round(length(recent_values)/2))-1])) / mean(recent_values) < 0.001
                elseif abs(maximum(recent_values[Int(round(length(recent_values) / 2)):end]) - maximum(recent_values[1:Int(round(length(recent_values) / 2))-1])) < 0.0001 &&
                    abs(minimum(recent_values[Int(round(length(recent_values) / 2)):end]) - minimum(recent_values[1:Int(round(length(recent_values) / 2))-1])) < 0.0001
                    println(2)
                    return true 
                end 
            end
            return false
        end

        # Affecter la condition d'arrêt
        function affect!(integrator)
            terminate!(integrator)
        end

        # Définir un callback pour arrêter l'intégration
        cb = DiscreteCallback(condition, affect!, save_positions=(false, false))

        #sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=frq_save,abstol=1e-9, reltol=1e-9, callback=cb) 
        sol = solve(prob, MethodOfSteps(RK4()), dt=pdt, saveat=frq_save) 

        Ò
        # Extraire les solutions
        t_values = sol.t # Les temps où la solution est évaluée
        
        u_values = sol.u # Les valeurs des solutions correspondantes
    elseif alg == "hard"
        # Définir le problème et le résoudre
        prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
        sol = solve(prob, MethodOfSteps(Vern9()), abstol=1e-9, reltol=1e-9, dt=1e-5)

        # Extraire les solutions
        t_values = sol.t # Les temps où la solution est évaluée
        u_values = sol.u # Les valeurs des solutions correspondantes
    else 
        # Définir le problème et le résoudre
        prob = DDEProblem(modele, u0, h, tspan, p; constant_lags = lags)
        sol = solve(prob, MethodOfSteps(Tsit5()))  
     
        # Extraire les solutions
        t_values = sol.t # Les temps où la solution est évaluée
        u_values = sol.u # Les valeurs des solutions correspondantes
    end