#####################################################
#### Analyse Modèle Ravageur 4 stades (Modèle 0) ####
#####################################################

# Charger les packages et les modules 

if Sys.KERNEL == :NT
    include("P:/Article1/fonction/fc_plots.jl")
    include("P:/Article1/fonction/fc_dynamique.jl")
else
    include("/Users/Quentin/Desktop/Article1/fonction/fc_plots.jl")
    include("/Users/Quentin/Desktop/Article1/fonction/fc_dynamique.jl")
end

using .module_graphique
using .module_analyse
using DifferentialEquations
using Plots
using NumericalIntegration
using Trapz
using Statistics
using DataFrames, CSV
using Dates
using Statistics


#######################################################
### Créaction des HeatMap Combiné pour les 7modèles ###
#######################################################

plotplot = false
localisation = "MeSUteste"
tstud = 500
iso = false


####################################################################
## Créaction des HeatMap Combiné pour le modèle sans parasitoïdes ##
####################################################################

module_graphique.heatmap_compose(iso, "MeSU", "sansPara", "DD", "q", "deltaU", tstud, "sansPara/LimiteStabDeltaU2")
module_analyse.printTemp()

module_graphique.heatmap_compose(iso, "MeSU", "sansPara", "DD", "q", "deltaL", tstud, "/sansPara/LimiteStabDeltaL2")
module_analyse.printTemp()

module_graphique.heatmap_compose(iso, "MeSU", "sansPara", "DD", "q", "deltaM", tstud, "/sansPara/LimiteStabDeltaM2")
module_analyse.printTemp()

module_graphique.heatmap_compose(iso, "MeSU", "sansPara", "DD", "q", "deltaA", tstud, "/sansPara/LimiteStabDeltaA1")
module_analyse.printTemp()

###########################################################################
## Créaction des HeatMap Combiné pour le modèle avec parasitoïde d'oeufs ##
###########################################################################

### Sans densité dépendence ###
module_graphique.heatmap_compose(iso, localisation, "ParaOeufs", "ssDD", "q", "deltaU", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaOeufs", "ssDD", "q", "deltaL", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaOeufs", "ssDD", "q", "deltaM", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaOeufs", "ssDD", "q", "deltaA", tstud)
module_analyse.printTemp(plotplot)



### Sans densité dépendence ###
module_graphique.heatmap_compose(iso, localisation, "ParaOeufs", "DD", "q", "deltaU", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaOeufs", "DD", "q", "deltaL", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaOeufs", "DD", "q", "deltaM", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaOeufs", "DD", "q", "deltaA", tstud)
module_analyse.printTemp(plotplot)

############################################################################
## Créaction des HeatMap Combiné pour le modèle avec parasitoïde de larve ##
############################################################################

### Sans densité dépendence ###
module_graphique.heatmap_compose(iso, localisation, "ParaLarve", "ssDD", "q", "deltaU", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaLarve", "ssDD", "q", "deltaL", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaLarve", "ssDD", "q", "deltaM", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaLarve", "ssDD", "q", "deltaA", tstud)
module_analyse.printTemp(plotplot)



### Sans densité dépendence ###
module_graphique.heatmap_compose(iso, localisation, "ParaLarve", "DD", "q", "deltaU", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaLarve", "DD", "q", "deltaL", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaLarve", "DD", "q", "deltaM", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaLarve", "DD", "q", "deltaA", tstud)
module_analyse.printTemp(plotplot)


#############################################################################
## Créaction des HeatMap Combiné pour le modèle avec parasitoïde de nymphe ##
#############################################################################

### Sans densité dépendence ###
module_graphique.heatmap_compose(iso, localisation, "ParaNymphe", "ssDD", "q", "deltaU", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaNymphe", "ssDD", "q", "deltaL", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaNymphe", "ssDD", "q", "deltaM", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaNymphe", "ssDD", "q", "deltaA", tstud)
module_analyse.printTemp(plotplot)



### Sans densité dépendence ###
module_graphique.heatmap_compose(iso, localisation, "ParaNymphe", "DD", "q", "deltaU", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaNymphe", "DD", "q", "deltaL", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaNymphe", "DD", "q", "deltaM", tstud)
module_analyse.printTemp(plotplot)

module_graphique.heatmap_compose(iso, localisation, "ParaNymphe", "DD", "q", "deltaA", tstud)
module_analyse.printTemp(plotplot)

