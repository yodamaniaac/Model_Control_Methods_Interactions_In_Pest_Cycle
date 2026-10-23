### Load necessary libraries ----
library(ggplot2)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(patchwork)
library(grid)
library(gridExtra)
library(cowplot)
library(ggsci)
library(scales)
library(colorspace)
library(svglite)
library(tidyr)


#################################
########### Fonctions ###########
#################################


### Base function for the creation of a single heatmap ----

# Parameters
#   pop_obs : character string
#       Prefix indicating the type of data to load
#       (e.g. "_heatmap_Densité_", "_heatmap_Variance_").
#
#   model_name : character string
#       Name of the model used to generate the data.
#
#   variable_name : character vector of length 2
#       Names of the parameters used as heatmap axes:
#           variable_name[1] : parameter on the Y-axis
#           variable_name[2] : parameter on the X-axis.
#
#   ligne_name : character string
#       Identifier of the parameter set (used for table organisation,
#       not directly in the plot).
#
#   variable_legende : character string
#       Biological compartment associated with the Y-axis
#       (e.g. egg, larva, pupa, adult).
#
#   isDD : character string
#       Flag indicating whether density dependence is included
#       (e.g. "DD" or "ssDD").
#
#   racine : character string
#       Root directory where CSV files are stored.
#
#   global_min, global_max : numeric
#       Global minimum and maximum values used to standardize
#       the color scale across multiple heatmaps.
#
#   Xtaille : numeric
#       Aspect ratio of the heatmap (passed to coord_fixed).
#
#   x_grille, y_grille : integer
#       Position of the plot within the global grid layout,
#       used to control axis labels and titles.
#
#   title_names : character string
#       Title displayed at the top of the heatmap when relevant.
#
#   plot_processing : character string, default = "identity"
#       Transformation applied to the color scale
#       (e.g. "identity", "log10").
#
#   plot_color : character string, default = "viridis"
#       Color palette used for the heatmap.
#
# Description
#   This function loads a precomputed CSV heatmap corresponding to a given
#   model and parameter pair, reshapes the data, and produces a ggplot2
#   heatmap with a standardized color scale. It is designed to be used
#   both standalone and as part of multi-panel comparison figures.
#
# Returns
#   A ggplot object representing the heatmap.
create_plot <- function(pop_obs ,model_name, variable_name, ligne_name, variable_legende, isDD, racine, global_min, global_max, 
                        Xtaille, x_grille, y_grille, title_names, plot_processing="identity", plot_color="viridis") {
  
  # Load the data
  file_path <- paste0(racine, model_name, pop_obs, isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  heatmap_densite <- read.csv(file_path)
  colnames(heatmap_densite) <- sub("^X", "", colnames(heatmap_densite))
  
  # Prepare heatmap data
  heatmap_melted <- melt(heatmap_densite, id.vars = "Row.Labels")
  colnames(heatmap_melted) <- c("Y", "X", "Density")
  colnames(heatmap_melted) <- c("Y", "X", "Density")
  heatmap_melted$X <- as.numeric(as.character(heatmap_melted$X))
  heatmap_melted$Density = round(heatmap_melted$Density, 7)
  heatmap_melted$X = round(heatmap_melted$X, 7)
  heatmap_melted$Y = round(heatmap_melted$Y, 7)
  
  
  
  # Fit the Y legend for additional mortality
  if (variable_name[1]=="deltaA")  {
    heatmap_melted$Y = heatmap_melted$Y - 0.1
  }
  
  
  # Create the base heatmap plot
  axeY <- paste0(variable_name[1], " : mortality of ",variable_legende)
  
  
  # Add 1 for the log10(value+1)
  if (plot_processing == "log10"){
    heatmap_melted$Density = heatmap_melted$Density + 1
  }
  
  # Set the name of the scale title
  if (pop_obs == "_heatmap_Variance_") {
    scale_title = "Variance"
  } else {
    scale_title = "Pest larval ensity"
  }
  
  #Creat the plot
  p <- ggplot(heatmap_melted) +
    geom_tile(aes(x = X, y = Y, fill = Density)) +
    scale_fill_viridis_c(name = scale_title ,trans = plot_processing, option = plot_color, limits = c(global_min,global_max)) +
    theme_minimal() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed(ratio = Xtaille) +
    labs( 
      title = if (y_grille == 1) paste(model_name, isDD) else NULL,
      x = if (y_grille == 4) "q : fecundity rate " else NULL,
      y = if (x_grille %% 4 == 1) axeY else NULL) +
    theme(
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  return(p)
}

### Base function for the creation of system state heatmaps ----

# Parameters
#   pop_obs : character string
#       Prefix indicating the type of population data to load
#       (typically pest density heatmaps).
#
#   model_name : character string
#       Name of the model used to generate the simulation outputs.
#
#   variable_name : character vector of length 2
#       Names of the parameters used as heatmap axes:
#           variable_name[1] : mortality parameter on the Y-axis
#           variable_name[2] : fecundity parameter on the X-axis.
#
#   ligne_name : character string
#       Identifier of the parameter set corresponding to the current row
#       in comparison tables (not directly used in the plot).
#
#   variable_legende : character string
#       Biological compartment associated with the mortality parameter
#       (e.g. egg, larva, pupa, adult).
#
#   isDD : character string
#       Flag indicating whether density dependence is included
#       (e.g. "DD" or "ssDD").
#
#   racine : character string
#       Root directory where CSV heatmap files are stored.
#
#   global_min, global_max : numeric
#       Global minimum and maximum values (not used directly here,
#       kept for consistency with other heatmap functions).
#
#   Xtaille : numeric
#       Aspect ratio of the heatmap (passed to coord_fixed).
#
#   x_grille, y_grille : integer
#       Position of the plot within the global grid layout,
#       used to control axis labels and titles.
#
#   title_name : character string
#       Title displayed at the top of the heatmap when relevant.
#
# Description
#   This function classifies the long-term dynamical state of the system
#   for each parameter combination into discrete categories:
#       - extinction,
#       - stable or unstable dynamics,
#       - presence or absence of parasitoids.
#   The classification is based on pest density, parasitoid density,
#   and variance thresholds computed from preprocessed CSV files.
#   The resulting heatmap provides a qualitative overview of system
#   stability across the parameter space.
#
# Returns
#   A ggplot object representing a categorical heatmap of system states,
#   including the pest survival boundary as a dashed isocline.
# Base function for the creation HeatMap states
create_plot_States <- function(pop_obs ,model_name, variable_name, ligne_name, variable_legende, isDD, racine, Xtaille, x_grille, y_grille, title_name) {
  
  # Find path and prepare pest larval density data
  file_path <- paste0(racine, model_name, pop_obs, isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  heatmap_densite <- read.csv(file_path)
  colnames(heatmap_densite) <- sub("^X", "", colnames(heatmap_densite))  
  heatmap_melted <- melt(heatmap_densite, id.vars = "Row.Labels")
  colnames(heatmap_melted) <- c("Y", "X", "Density")
  heatmap_melted$X <- as.numeric(as.character(heatmap_melted$X))
  
  
  # Find path and prepare variance data 
  file_path <- paste0(racine, model_name, '_heatmap_Variance_',isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  heatmap_Variance <- read.csv(file_path)
  colnames(heatmap_Variance) <- sub("^X", "", colnames(heatmap_Variance))
  isocline_melted <- melt(heatmap_Variance, id.vars = "Row.Labels")
  colnames(isocline_melted) <- c("Y", "X", "Density")
  isocline_melted$X <- as.numeric(as.character(isocline_melted$X))
  
  
  # Find path and prepare adult parasitoid density data
  file_path <- paste0(racine, model_name, '_heatmap_Densité_Para_',isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  heatmap_densite_para <- read.csv(file_path)
  colnames(heatmap_densite_para) <- sub("^X", "", colnames(heatmap_densite_para))
  heatmap_para_melted <- melt(heatmap_densite_para, id.vars = "Row.Labels")
  colnames(heatmap_para_melted) <- c("Y", "X", "Density")
  heatmap_para_melted$X <- as.numeric(as.character(heatmap_para_melted$X))
  
  # Fit the Y legend for additional mortality
  if (variable_name[1]=="deltaA")  {
    heatmap_melted$Y = heatmap_melted$Y - 0.1
    heatmap_para_melted$Y = heatmap_para_melted$Y - 0.1
    isocline_melted$Y = isocline_melted$Y - 0.1
  }
  
  
  
  # Create existence limit
  taille = 300
  X = seq(from = min(isocline_melted$X), to = max(isocline_melted$X), length.out = taille)
  Y = rep(0,taille)
  inc = 0
  for (i in X) {
    inc = inc + 1 
    if (variable_name[1]=="deltaU"){
      Y[inc] = -(1/4) * log(0.1/i)
    }else if (variable_name[1]=="deltaL"){
      Y[inc] = -(1/16) * log(0.1/i)   
    }else if (variable_name[1]=="deltaM"){
      Y[inc] = -(1/8) * log(0.1/i)
    }
  }
  isocline_existance <- data.frame("X" = X, "Y" = Y)
  isocline_existance <- isocline_existance[isocline_existance$Y <= 1, ]
  
  
  
  # Name Y axe
  axeY <- paste0(" mortality rate of ",variable_legende)
  
  # Create system status categories
  parasitoid_present <- heatmap_para_melted$Density > 1e-6
  host_alive <- heatmap_melted$Density > 1
  host_stable <- isocline_melted$Density < 0.01
  
  for (i in seq_len(nrow(heatmap_para_melted))) {
    
    if (!host_alive[i]) {
      heatmap_para_melted$Density[i] <- "All dead"
      
    } else if (!parasitoid_present[i]) {
      
      if (host_stable[i]) {
        heatmap_para_melted$Density[i] <- "Stable without parasitoid"
      } else {
        heatmap_para_melted$Density[i] <- "Unstable without parasitoid"
      }
      
    } else {
      
      if (host_stable[i]) {
        heatmap_para_melted$Density[i] <- "Stable with parasitoid"
      } else {
        heatmap_para_melted$Density[i] <- "Unstable with parasitoid"
      }
    }
  }
  
  heatmap_para_melted$Density <- factor(heatmap_para_melted$Density, levels = c("All dead", "Stable without parasitoid", "Unstable without parasitoid", "Stable with parasitoid", "Unstable with parasitoid"))
  
  
  ## Create the plot 
  p <- ggplot(heatmap_para_melted) +
    geom_tile(aes(x = X, y = Y, fill = Density)) +
    scale_fill_manual(
      name = "System state",
      
      # Choose the color
      values = c(
        "All dead"= "#E3F2FD",
        "Stable without parasitoid" =  "#B5C861FF",
        "Unstable without parasitoid" = "#8AA789FF",
        "Stable with parasitoid" =  "#CB8573FF",
        "Unstable with parasitoid" = "#976153FF"
      )
    ) +
    theme_minimal() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed(ratio = Xtaille) +
    labs( 
      title = if (y_grille == 1) paste(title_name) else NULL,
      x = if (y_grille == 4) "fecondity rate" else NULL,
      y = if (x_grille %% 4 == 1) axeY else NULL
    ) +
    theme(
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  # Add the area on pest survival
  if (variable_name[1]!="deltaA") {
    p = p + geom_line(data = isocline_existance, aes(x = X, y = Y), 
                      color = "#666666", linetype = "dashed", linewidth = 1)
  }
  
  return(p)
}


### Base function for the creation of relative effect heatmaps ----

# Parameters
#   pop_obs : character string
#       Prefix indicating the type of population data to load
#       (e.g. pest density heatmaps).
#
#   model_name : character string
#       Name of the model including parasitoids.
#
#   variable_name : character vector of length 2
#       Names of the parameters used as heatmap axes:
#           variable_name[1] : mortality parameter on the Y-axis
#           variable_name[2] : fecundity parameter on the X-axis.
#
#   ligne_name : character string
#       Identifier of the parameter set corresponding to the current row
#       in comparison tables (not directly used in the plot).
#
#   variable_legende : character string
#       Biological compartment associated with the mortality parameter
#       (e.g. egg, larva, pupa, adult).
#
#   isDD : character string
#       Flag indicating whether density dependence is included
#       (e.g. "DD" or "ssDD").
#
#   racine : character string
#       Root directory where CSV heatmap files are stored.
#
#   global_min, global_max : numeric
#       Global minimum and maximum values (currently not used,
#       kept for consistency with other heatmap functions).
#
#   Xtaille : numeric
#       Aspect ratio of the heatmap (passed to coord_fixed).
#
#   x_grille, y_grille : integer
#       Position of the plot within the global grid layout,
#       used to control axis labels and titles.
#
#   title_names : character string
#       Title displayed at the top of the heatmap when relevant.
#
# Description
#   This function computes and visualizes the relative effect of parasitoids
#   on pest population density by comparing simulations with and without
#   parasitoids. The heatmap represents the difference in pest density
#   between the two scenarios across the parameter space.
#   Positive values indicate a reduction in pest density due to parasitoids.
#
# Returns
#   A ggplot object representing a continuous heatmap of relative effects.
# Base function for the creation HeatMap states
create_plot_Relatif <- function(pop_obs ,model_name, variable_name, ligne_name, variable_legende, isDD, racine, global_min, global_max, Xtaille, x_grille, y_grille, title_names) {
  
  # Alternative Para
  file_path <- paste0(racine, model_name, pop_obs, isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  heatmap_densite <- read.csv(file_path)
  colnames(heatmap_densite) <- sub("^X", "", colnames(heatmap_densite))
  
  file_path <- paste0(racine, model_name, '_heatmap_Variance_',isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  heatmap_Variance <- read.csv(file_path)
  colnames(heatmap_Variance) <- sub("^X", "", colnames(heatmap_Variance))
  
  # Prepare with out parasitoid data
  file_path <- paste0(racine, "sansPara", pop_obs, isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  if (file.exists(file_path)) {
    heatmap_densite_without_para <- read.csv(file_path)
    colnames(heatmap_densite_without_para) <- sub("^X", "", colnames(heatmap_densite_without_para))
    heatmap_densite_without_para_melted <- melt(heatmap_densite_without_para, id.vars = "Row.Labels")
    colnames(heatmap_densite_without_para_melted) <- c("Y", "X", "Density")
    heatmap_densite_without_para_melted$X <- as.numeric(as.character(heatmap_densite_without_para_melted$X))
  }
  
  
  # Prepare heatmap data
  heatmap_melted <- melt(heatmap_densite, id.vars = "Row.Labels")
  colnames(heatmap_melted) <- c("Y", "X", "Density")
  heatmap_melted$X <- as.numeric(as.character(heatmap_melted$X))
  
  # Prepare isocline data
  isocline_melted <- melt(heatmap_Variance, id.vars = "Row.Labels")
  colnames(isocline_melted) <- c("Y", "X", "Density")
  isocline_melted$X <- as.numeric(as.character(isocline_melted$X))
  
  
  # Create the base heatmap plot
  axeY <- paste0(variable_name[1], " : mortalité des ",variable_legende)
  
  #Preparer les données transformé
  heatmap_plot <- data.frame(
    X = heatmap_melted$X,
    Y = heatmap_melted$Y,
    #Density = pmax(heatmap_melted$Density / heatmap_densite_without_para_melted$Density, 0)
    Density = heatmap_densite_without_para_melted$Density - heatmap_melted$Density 
  )
  
  # Ta HeatMap avec couleurs définies manuellement
  p <- ggplot(heatmap_plot) +
    geom_tile(aes(x = X, y = Y, fill = Density)) +
    #scale_fill_viridis_c(limits = c(global_min, global_max), 
    #      guide = if (y_grille == 1) "colorbar" else "none") +
    scale_fill_viridis_c(limits = c(0, 1)) +
    #scale_fill_viridis_c() +
    #scale_fill_viridis_c(trans = "log10") +
    theme_minimal() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed(ratio = Xtaille) +
    labs( 
      title = if (y_grille == 1) paste(model_name, isDD) else NULL,
      x = if (y_grille == 4) "q : taux de reproduction" else NULL,
      y = if (x_grille %% 4 == 1) axeY else NULL) +
    theme(
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  return(p)
}




### Utility function to compute global color scale limits ----

# Parameters
#   pop_obs : character string
#       Prefix indicating the type of heatmap data to load
#       (e.g. "_heatmap_Densité_", "_heatmap_Variance_").
#
#   racine : character string
#       Root directory where CSV heatmap files are stored.
#
#   model_names : character vector
#       Names of the models for which heatmap data are available.
#
#   variable_names1 : character vector
#       Names of the parameters used on the Y-axis of the heatmaps.
#
#   variable_names2 : character vector
#       Names of the parameters used on the X-axis of the heatmaps.
#
#   isDD : character string
#       Flag indicating whether density dependence is included
#       (e.g. "DD" or "ssDD").
#
# Description
#   This function scans all available heatmap CSV files corresponding
#   to the provided models and parameter combinations in order to find
#   the global minimum and maximum values. These values are typically
#   used to enforce a consistent color scale across multiple heatmaps
#   and comparison figures.
#
# Returns
#   A named list with two elements:
#       - min : global minimum value across all heatmaps
#       - max : global maximum value across all heatmaps
find_global_min_max <- function(pop_obs ,racine, model_names, variable_names1, variable_names2, isDD) {
  
  # Set maximum and minimum value 
  min_val <- Inf
  max_val <- -Inf
  
  # Look all the value to find the min and max
  for (model_name in model_names) {
    for (variable_name1 in variable_names1) {
      for (variable_name2 in variable_names2) {
        # Path to the File
        file_path <- paste0(racine, model_name, pop_obs, isDD,"_", variable_name2,"_", variable_name1, '.csv')
        
        # Look the values if the file exist
        if (file.exists(file_path)) {
          heatmap_densite <- read.csv(file_path)
          heatmap_densite <- heatmap_densite[, -1]  # Remove row labels
          heatmap_densite <- as.matrix(heatmap_densite)
          heatmap_densite <- as.numeric(heatmap_densite)
          
          current_min <- min(heatmap_densite, na.rm = TRUE)
          current_max <- max(heatmap_densite, na.rm = TRUE)
          
          # Update the min and max
          if (current_min < min_val) min_val <- current_min
          if (current_max > max_val) max_val <- current_max
          
          
        }
      }
    }
  }
  
  return(list(min = min_val, max = max_val))
}




### Main function to generate comparison tables of heatmaps ----

# Parameters
#   Empfile : character string
#       Name of the empirical dataset or experiment used to build
#       file paths and output directories.
#
#   type : character string, default = "Density"
#       Type of quantity visualised in the heatmaps.
#       Supported values include:
#           - "Density" : population density
#           - "CV"      : coefficient of variation
#           - "Type"    : dominant state / regime
#           - "Relatif" : relative values between scenarios
#
#   variable_names : character vector
#       Names of the primary parameters varied along the Y-axis
#       of the heatmaps (e.g. mortality rates for each stage).
#
#   variable2_names : character string
#       Name of the secondary parameter varied along the X-axis
#       of the heatmaps (typically "q").
#
#   ligne_names : character vector
#       Labels used to identify each row of the comparison table,
#       corresponding to variable_names.
#
#   variable_legende : character vector
#       Names of the biological compartments displayed in the legend
#       (e.g. egg, larva, pupa, adult).
#
#   title_names : character vector
#       Titles of the columns, corresponding to the different
#       parasitoid scenarios.
#
#   model_names : character vector
#       Internal model identifiers used to locate the heatmap files
#       (must match file naming conventions).
#
# Description
#   This function builds a full comparison table of heatmaps across
#   multiple models and parameters under density-dependent (DD)
#   conditions. It ensures a consistent color scale across all panels
#   by computing global minimum and maximum values beforehand.
#
#   The resulting figure is organised as:
#       - columns: parasitoid scenarios
#       - rows   : life-stage-specific parameters
#       - shared legend and common axes
#
#   The function automatically adapts file paths depending on the
#   operating system and saves the final figure in both SVG and PNG
#   formats.
#
# Returns
#   A patchwork / ggplot object corresponding to the full comparison
#   table, invisibly printed and saved to disk.
create_comparison_tables <- function(Empfile, type = "Density",variable_names = c("deltaU","deltaL","deltaM","deltaA"),variable2_names="q",ligne_names = c("DeltaA1","DeltaU1","DeltaL2","DeltaM2"),
                                     variable_legende = c("egg", "larva", "pupa", "adult"), title_names = c("Without parasitoid","Egg parasitoid","Larva parasitoid","Pupa parasitoid"), 
                                     model_names = c("sansPara", "ParaOeufs", "ParaLarve", "ParaNymphe")) {
  
  
  # Check the operating system and set the root directory accordingly
  if (Sys.info()["sysname"] == "Windows") {
    racine <-  paste0("P:/Article1/tables/",Empfile,"/")
    racineSave <-  paste0("P:/Article1/tables/graphs/graphR/",Empfile,"/")
  } else {
    racine <- paste0("/Users/Quentin/Desktop/Article1/tables/",Empfile,"/")
    racineSave <- paste0("/Users/Quentin/Desktop/Article1/graphs/graphR/",Empfile,"/")
  }
  
  # We only plot the Pest larval density
  if (type == "CV"){
    pop_obs <- "_heatmap_Variance_"
  } else {
    pop_obs <- "_heatmap_Densité_"   
  }
  
  
  
  
  # Find global min and max for DD and ssDD
  global_min_max_DD <- find_global_min_max(pop_obs, racine, model_names, variable_names[1:4], variable2_names,"DD")
  Plot <- list()
  
  
  o = 0
  count_col = 1
  for (i in 1:length(variable_names)) {
    variable_name = variable_names[i]
    ligne_name = ligne_names[i]
    o = o+1
    model_name = model_names[1]
    if (type == "Density" && variable2_names == "q"){
      
      Plot[[paste0(model_name, "_", variable_name)]] <- create_plot(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                    "DD", racine, global_min_max_DD$min, global_min_max_DD$max, 70, o, i, title_names[count_col],"log10","viridis")
      
    } else if (type == "Type" && variable2_names == "q") {
      Plot[[paste0(model_name, "_", variable_name)]] <- create_plot_States(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                           "DD", racine, 70, o, i, title_names[count_col])
      
    } else if (type == "CV" && variable2_names == "q") {
      Plot[[paste0(model_name, "_", variable_name)]] <- create_plot(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                    "DD", racine, global_min_max_DD$min, global_min_max_DD$max, 70, o, i, title_names[count_col],"identity","cividis")
    }
    
    for (model_name in model_names[2:length(model_names)]) {
      o = o+1
      count_col = count_col + 1
      
      if (type == "Density"){
        # Generate plots for DD
        Plot[[paste0(model_name, "_", variable_name)]] <- create_plot(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                      "DD", racine, global_min_max_DD$min, global_min_max_DD$max, 70, o, i, title_names[count_col],"log10","viridis")
      } else if (type == "Type") {
        Plot[[paste0(model_name, "_", variable_name)]] <- create_plot_States(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                             "DD", racine, 70, o, i, title_names[count_col])
      } else if (type == "CV"&& variable2_names == "q") {
        Plot[[paste0(model_name, "_", variable_name)]] <- create_plot(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                      "DD", racine, global_min_max_DD$min, global_min_max_DD$max, 70, o, i, title_names[count_col],"identity","cividis")
      }else if (type == "Relatif") {
        Plot[[paste0(model_name, "_", variable_name)]] <- create_plot_Relatif(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                              "DD", racine, global_min_max_DD$min, global_min_max_DD$max, 70, o, i, title_names[count_col])
      }
    }
  }
  
  
  
  
  # Récupérer la légende depuis un des plots d'origine
  legend_plot <- Plot[[10]] + 
    theme(
      legend.position = "right",
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 30, face = "bold"),
      legend.key.size = unit(1.2, "cm"),
      legend.spacing.y = unit(0.5, "cm"),
      panel.background  = element_rect(fill = hcl(39, 100, 95), colour = NA),
      plot.background   = element_rect(fill = hcl(39, 100, 95), colour = NA),
    )
  legend <- cowplot::get_legend(legend_plot )
  
  # Convertir en un ggdraw pour l'utiliser dans patchwork
  legend_gg <- cowplot::ggdraw(legend)
  
  
  # Création d'un objet de même dimension que la légende mais vide
  empty_legend <- Plot[[10]] + 
    theme(
      legend.position = "right",
      legend.text = element_text(size = 14, color =  "white"),
      legend.title = element_text(size = 16, face = "bold", color =  "white"),
      legend.key = element_rect(fill =  "white", colour =  "white"),  # fond des carrés de légende
      legend.background = element_rect(fill =  "white", colour =  "white"),  # fond du bloc de légende
      legend.box.background = element_rect(fill =  "white", colour =  "white"),  # bordure autour de toute la légende
      legend.key.size = unit(0, "cm"),
      legend.spacing.y = unit(0, "cm")
    )
  empty_legend <- cowplot::get_legend(empty_legend )
  empty_legend_gg <- cowplot::ggdraw(empty_legend)
  
  # Nettoyer tous les plots
  Plot_clean <- lapply(Plot, function(p) {
    p + 
      theme(
        legend.position = "none",
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      )
  })
  
  
  # Créer titres colonnes et lignes
  col_titles <- c("Without parasitoid", "Egg parasitoid", "Larva parasitoid", "Pupa parasitoid")
  row_labels <- c("Egg mortality rate", "Larva mortality rate", "Pupa mortality rate", "Adult mortality rate")
  
  # Créer chaque rangée de graphiques avec titre de ligne
  plot_rows <- list()
  for (i in 0:3) {
    label_grob <- ggplot() +
      theme_void() +
      annotation_custom(
        grob = textGrob(row_labels[i + 1], rot = 90, gp = gpar(fontsize = 30)),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      )
    row_plots <- Plot_clean[(1 + i*4):(4 + i*4)]
    
    # La légende uniquement pour la 1ère ligne
    legend_gg <- if (i == 0) legend_gg else empty_legend_gg
    
    row <- label_grob + row_plots[[1]] + row_plots[[2]] + row_plots[[3]] + row_plots[[4]] + legend_gg
    plot_rows[[i + 1]] <- row + plot_layout(ncol = 6, widths = c(0.1, rep(1, 4), 1.1)) # 7 colonne, 4 de 1 de large (les schéma), une de 0.1 et une de 0.2 (la légende)
  }
  
  # Ligne du haut avec titres de colonnes
  col_title_plots <- lapply(col_titles, function(title) {
    ggplot() +
      theme_void() +
      annotation_custom(textGrob(title, gp = gpar(fontsize = 30, fontface = "bold")))
  })
  empty_plot <- ggplot() + theme_void()
  empty_legend_spot <- ggplot() + theme_void()  # case vide pour la légende
  top_row <- empty_plot + col_title_plots[[1]] + col_title_plots[[2]] + col_title_plots[[3]] + col_title_plots[[4]] + empty_legend_spot 
  top_row <- top_row + plot_layout(ncol = 6, widths = c(0.1, rep(1, 4), 1.1))
  
  # Axe X commun
  bottom_axis <- ggplot() +
    theme_void() +
    annotation_custom(
      grob = textGrob("Fecondity rate", gp = gpar(fontsize = 30))
    )
  
  # Fusion finale
  final_plot <- top_row /
    plot_spacer() /
    plot_rows[[1]] /
    plot_rows[[2]] /
    plot_rows[[3]] /
    plot_rows[[4]] /
    bottom_axis
  
  final_plot <- final_plot +
    plot_layout(heights = c(0.1, 0.01, rep(1, 4), 0.1))
  
  # Affichage
  print(final_plot)
  
  # Sauvegarde
  
  plot_file_path_DD_svg <- paste0(racineSave, pop_obs, "_", type, "_", "comparison_table_DD.svg")
  plot_file_path_DD_png <- paste0(racineSave, pop_obs, "_", type, "_", "comparison_table_DD.png")
  
  ggsave(plot_file_path_DD_svg, device = "svg", plot = final_plot, width = 25, height = 23, dpi = 300,   bg ="#FFF6E5")
  ggsave(plot_file_path_DD_png, device = "png", plot = final_plot, width = 25, height = 23, dpi = 300,   bg ="#FFF6E5")
  
  cat("Comparison tables saved as", plot_file_path_DD_png, "and", plot_file_path_DD_svg, "\n")
  
  return(final_plot)
  
}











### Start creating figures  ----


## Plot for article 

# HeatMap density Final (dev des para 16) 
create_comparison_tables("HeatMapDD_paradev_16_final")
create_comparison_tables("Heatmap_para_16_morta_add")
create_comparison_tables("Heatmap_para_16_ss_competition")

### HeatMap densité Relative Final (dev des para 16)
create_comparison_tables("HeatMapDD_paradev_16_final", type = "Relatif")
create_comparison_tables("Heatmap_para_16_morta_add", type = "Relatif")

### HeatMap Type Final (dev des para 16) version type d'équilibre
create_comparison_tables("HeatMapDD_paradev_16_final", type = "Type")
create_comparison_tables("Heatmap_para_16_morta_add", type = "Type")

### HeatMap Coéficcient de variation Final (dev des para 16) version type d'équilibre (heuuu je crois que c'est plutôt la variance )
create_comparison_tables("HeatMapDD_paradev_16_final", type = "CV")
create_comparison_tables("Heatmap_para_16_morta_add", type = "CV")

### HeatMap densité Final (dev des para 16) version type d'équilibre
create_comparison_tables("HeatMapDD_paradev_DevMorta",FALSE, variable_names = c("deltaU","deltaL","deltaM","deltaA"),variable2_names="tauJ",
                         ligne_names = c("DeltaA1","DeltaU1","DeltaL2","DeltaM2"), variable_legende = c("oeufs", "larves", "nymphes","adultes"))
create_comparison_tables("Heatmap_para_16_morta_add",FALSE, variable_names = c("deltaU","deltaL","deltaM","deltaA"),variable2_names="tauJ",
                         ligne_names = c("DeltaA1","DeltaU1","DeltaL2","DeltaM2"), variable_legende = c("oeufs", "larves", "nymphes","adultes"))










### The Great Verification ----


#### attack rate

###### a = 0
create_comparison_tables("GrandeVerification/Verif_heatmap_a0/fixeSimple/0.0001",FALSE, "Type")
###### a = 0.03
create_comparison_tables("GrandeVerification/Verif_heatmap_a003/fixeSimple/0.0001",FALSE, "Type")
###### a = 0.001
create_comparison_tables("GrandeVerification/Verif_heatmap_a0001/fixeSimple/0.0001",FALSE, "Type")



#### Live time of pest juvenile stages = 4
create_comparison_tables("GrandeVerification/Verif_heatmap_tpsdev4/fixeSimple/0.0001",FALSE, "Type")



#### Pest juvenils mortality rate : 0.1
create_comparison_tables("GrandeVerification/Verif_heatmap_delta01/fixeSimple/0.0001",FALSE, "Type")



#### Conpetition variation

###### alphaL = 0.00001
create_comparison_tables("GrandeVerification/Verif_heatmap_alphaL000001/fixeSimple/0.0001",FALSE, "Type")
###### alphaL = 0.0000001
create_comparison_tables("GrandeVerification/Verif_heatmap_alphaL00000001/fixeSimple/0.0001",FALSE, "Type")
###### alphaL = 0
create_comparison_tables("GrandeVerification/Heatmap_para_16_ss_competition/fixeSimple/0.0001",FALSE, "Type")


#### Competition exponent #don't calculate until the end, look at the LOGS

###### exponent = 1.1
create_comparison_tables("GrandeVerification/Verif_heatmap_DD11/fixeSimple/0.0001",FALSE, "Type")
###### exponent = 1.5
create_comparison_tables("GrandeVerification/Verif_heatmap_DD15/fixeSimple/0.0001",FALSE, "Type")
###### exponent = 1.8
create_comparison_tables("GrandeVerification/Verif_heatmap_DD18/fixeSimple/0.0001",FALSE, "Type")



#### Variation in pest development time

###### tauL = 1
create_comparison_tables("GrandeVerification/Verif_heatmap_tauJ1/fixeSimple/0.0001",FALSE, "Type")
###### tauL = 14
create_comparison_tables("GrandeVerification/Verif_heatmap_tauJ14/fixeSimple/0.0001",FALSE, "Type")
###### tauL = 28
create_comparison_tables("GrandeVerification/Verif_heatmap_tauJ28/fixeSimple/0.0001",FALSE, "Type")


#### Type II RF: variation in handling time

###### H = 0
create_comparison_tables("GrandeVerification/Verif_heatmap_RFII0/fixeSimple/0.0001",FALSE, "Type")
###### H = 0.1
create_comparison_tables("GrandeVerification/Verif_heatmap_RFII01/fixeSimple/0.0001",FALSE, "Type")
###### H = 0.25
create_comparison_tables("GrandeVerification/Verif_heatmap_RFII025/fixeSimple/0.0001",FALSE, "Type")
###### H = 0.5
create_comparison_tables("GrandeVerification/Verif_heatmap_RFII05/fixeSimple/0.0001",FALSE, "Type")
###### H = 1
create_comparison_tables("GrandeVerification/Verif_heatmap_RFII1/fixeSimple/0.0001",FALSE, "Type")
###### H = 2
create_comparison_tables("GrandeVerification/Verif_heatmap_RFII2/fixeSimple/0.0001",FALSE, "Type")



