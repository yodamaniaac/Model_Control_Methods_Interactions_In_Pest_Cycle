# Model description

This GitHub repository presents a model of delay differential equations that models a crop pest population and various control methods. 


The "Mathematica" folder contains a semi-analytical study used to determine the existence limit of the pest in the absence of larval competition and the stability limit of the system with larval competition but without parasitoids. 


The "Simulation_Heatmaps" folder contains Julia code for simulating the heatmaps presented in the article. Bash scripts and a Julia environment are provided to run these codes on a computer in a parallelized manner.
The results of these simulations are located in the "table" folder. 
This data is then plotted using the R code `fc_graph_tot.R` located in the “graphs” folder. All the graphs are also located in this same folder, in the “Final_Graph” directory.
The other graphs—the bifurcation diagram and the one-dimensional figure—are simulated and plotted using the Julia code found in “Local_Simulation.”

# Heatmaps simulation
Creating heatmaps involves a series of resource-intensive simulations.  These tasks are therefore parallelized using the following Bash script: bashHeatMapDD.sh  
This script calls Main_Modele0_MeSU.jl, Main_Modele1_DD_MeSU.jl, Main_Modele2_DD_MeSU.jl, Main_Modele3_DD_MeSU.jl, which launch the various simulations for models without parasitoids, with parasitoids, and for larvae, eggs, and nymphs.  
All of these files call the function file fc_dynamique.jl, which provides the various functions needed to simulate population dynamics.  
The simulated data are stored in the "table" folder.  
This data is used to plot various heatmaps representing the density of pest larvae and adult parasitoids, the different zones of stability and survival for the pest, and the coefficient of variation.  
These graphs are located in "graphs/Final_Graph/Heatmaps."
The graph showing the stability and survival zones of the parasitoid is plotted on top of the stability limits for the parasitoid and the survival limits for the pest, all of which were calculated in Mathematica files.
# Other graph simulation
Bifurcation diagrams and one-dimensional figures require less computational power and can therefore be easily calculated locally using the code provided in Local_Simulation.
fc_local_simulation provides all the functions necessary for simulating and plotting these figures for the various models studied.  
BifurcationFigH.jl allows for the plotting of bifurcation curves for the densities of pest larvae, adult pests, and adult parasitoids as a function of larval mortality rates for a parasitoid model involving nymphs.  
It is linked to BifurcationSurvieNymphes.jl, which plots nymph survival rates under the same conditions.  
Para1D.jl is used to plot whether the parasitoid survives or not based on the intensity of larval competition.  
All of these figures can be found in "graphs/Final_Graph/Bifurcation" or "Figure1D"