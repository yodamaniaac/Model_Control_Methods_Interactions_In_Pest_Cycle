# Load necessary libraries
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
#################################
########### Fonctions ###########
#################################


# Base function for the creation of a signle Heatmap
create_plot <- function(pop_obs ,model_name, variable_name, ligne_name, variable_legende, isDD, racine, isLimite, global_min, global_max, Xtaille, x_grille, y_grille, title_names) {
  # Load the data
  # Densité Oeufs
  #file_path <- paste0(racine, model_name, '_heatmap_Densité_',isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  
  # Alternative Para
  file_path <- paste0(racine, model_name, pop_obs, isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  
  
  heatmap_densite <- read.csv(file_path)
  colnames(heatmap_densite) <- sub("^X", "", colnames(heatmap_densite))
  
  file_path <- paste0(racine, model_name, '_heatmap_Variance_',isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  heatmap_Variance <- read.csv(file_path)
  colnames(heatmap_Variance) <- sub("^X", "", colnames(heatmap_Variance))
  
  if (isLimite == TRUE && isDD == "DD") {
    if (Sys.info()["sysname"] == "Windows") {
      racineLim <-  paste0("P:/Article1/tables/limite/")
    } else {
      racineLim <- paste0("/Users/Quentin/Desktop/Article1/tables/limite/")
    }
    file_path <- paste0(racineLim, model_name, '_LimiteStab',ligne_name,'.csv')
    heatmap_Limite <- read.csv(file_path, header = TRUE)
    colnames(heatmap_Limite) <- sub("^X", "", colnames(heatmap_Limite))
    heatmap_Limite <- rbind(colnames(heatmap_Limite), heatmap_Limite)
    colnames(heatmap_Limite) <- c("X", "Y")
    heatmap_Limite$X <- as.numeric(as.character(heatmap_Limite$X))
    heatmap_Limite$Y <- as.numeric(as.character(heatmap_Limite$Y))
    heatmap_Limite <- heatmap_Limite[heatmap_Limite$X <= Xtaille, ]
    heatmap_Limite$Density = pmax(heatmap_Limite$Density, 10^-10)
  }
  
  # Prepare heatmap data
  heatmap_melted <- melt(heatmap_densite, id.vars = "Row.Labels")
  colnames(heatmap_melted) <- c("Y", "X", "Density")
  heatmap_melted$X <- as.numeric(as.character(heatmap_melted$X))
  #heatmap_melted$Density = pmax(heatmap_melted$Density, 10^-10)
  heatmap_melted$Density = heatmap_melted$Density + 1
  
  # Prepare isocline data
  isocline_melted <- melt(heatmap_Variance, id.vars = "Row.Labels")
  colnames(isocline_melted) <- c("Y", "X", "Density")
  isocline_melted$X <- as.numeric(as.character(isocline_melted$X))
  isocline_melted$Density = pmax(isocline_melted$Density, 10^-10)
  
  
  # Create the base heatmap plot
  axeY <- paste0(variable_name[1], " : mortalité des ",variable_legende)
  if (isDD=="ssDD") {
    p <- ggplot(heatmap_melted) +
      geom_tile(aes(x = X, y = Y, fill = Density)) +
      scale_fill_viridis_c(limits = c(global_min, global_max)) +
      theme_minimal() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed(ratio = Xtaille) +
      labs( 
        title = if (y_grille == 1) paste(model_name, isDD) else NULL,
        x = if (y_grille == 4) "q : taux de reproduction" else NULL,
        y = if (x_grille %% 4 == 1) axeY else NULL)
  }else{
    p <- ggplot(heatmap_melted) +
      geom_tile(aes(x = X, y = Y, fill = Density)) +
      #scale_fill_viridis_c(limits = c(global_min, global_max), 
      #      guide = if (y_grille == 1) "colorbar" else "none") +
      scale_fill_viridis_c(trans = "log10", limits = c(global_min,global_max)) +
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
  }
  
  
  # Add isoclines
  #q <- p + geom_contour(data = isocline_melted, aes(x = X, y = Y, z = Density),
  #                      color ="yellow", alpha = 1, breaks = c(1), linewidth = 1) 
  #geom_contour(data = isocline_melted, aes(x = X, y = Y, z = Density),
  #             color ="orange", alpha = 1, breaks = c(100), linewidth = 1) +
  #geom_contour(data = isocline_melted, aes(x = X, y = Y, z = Density),
  #             color ="red", alpha = 1, breaks = c(200), linewidth = 1)+
  #geom_contour(data = heatmap_melted, aes(x = X, y = Y, z = Density),
  #             color ="purple", alpha = 1, breaks = c(0.0001), linewidth = 1)
  
  # Add isocline for para
  #file_path <- paste0(racine, model_name, '_heatmap_Densité_Para_',isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  #if (file.exists(file_path)) {
  #  heatmap_densite_para <- read.csv(file_path)
  #  colnames(heatmap_densite_para) <- sub("^X", "", colnames(heatmap_densite_para))
  #  heatmap_para_melted <- melt(heatmap_densite_para, id.vars = "Row.Labels")
  #  colnames(heatmap_para_melted) <- c("Y", "X", "Density")
  #  heatmap_para_melted$X <- as.numeric(as.character(heatmap_para_melted$X))  
  #  heatmap_para_melted$Density = pmax(heatmap_para_melted$Density, 10^-10)
  #  q <- q + geom_contour(data = heatmap_para_melted, aes(x = X, y = Y, z = Density),
  #                        color ="cyan", alpha = 1, breaks = c(0.0001), linewidth = 1)}
  
  
  #if (isLimite == TRUE) {
  # Add stability limit curve
  #  heatmap_Limite <- heatmap_Limite[order(heatmap_Limite$Y), ]
  #  heatmap_Limite$Density = pmax(heatmap_Limite$Density, 10^-10)
  #  r <- q + geom_point(data = heatmap_Limite, aes(y = Y, x = X), color = "grey", size = 0.01, alpha = 10)
  #} else {
  #  r <- q
  #}
  
  return(p)
}


# Base function for the creation HeatMap states
create_plot_States <- function(pop_obs ,model_name, variable_name, ligne_name, variable_legende, isDD, racine, isLimite, global_min, global_max, Xtaille, x_grille, y_grille, title_name) {
  
  # Load the data
  # Densité Oeufs
  #file_path <- paste0(racine, model_name, '_heatmap_Densité_',isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  
  # Alternative Para
  file_path <- paste0(racine, model_name, pop_obs, isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  
  
  heatmap_densite <- read.csv(file_path)
  colnames(heatmap_densite) <- sub("^X", "", colnames(heatmap_densite))
  
  file_path <- paste0(racine, model_name, '_heatmap_Variance_',isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  heatmap_Variance <- read.csv(file_path)
  colnames(heatmap_Variance) <- sub("^X", "", colnames(heatmap_Variance))
  
  # Add isocline for para
  file_path <- paste0(racine, model_name, '_heatmap_Densité_Para_',isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  if (file.exists(file_path)) {
    heatmap_densite_para <- read.csv(file_path)
    colnames(heatmap_densite_para) <- sub("^X", "", colnames(heatmap_densite_para))
    heatmap_para_melted <- melt(heatmap_densite_para, id.vars = "Row.Labels")
    colnames(heatmap_para_melted) <- c("Y", "X", "Density")
  }
  
  
  # Prepare heatmap data
  heatmap_melted <- melt(heatmap_densite, id.vars = "Row.Labels")
  colnames(heatmap_melted) <- c("Y", "X", "Density")
  heatmap_melted$X <- as.numeric(as.character(heatmap_melted$X))
  
  # Prepare isocline data
  isocline_melted <- melt(heatmap_Variance, id.vars = "Row.Labels")
  colnames(isocline_melted) <- c("Y", "X", "Density")
  isocline_melted$X <- as.numeric(as.character(isocline_melted$X))
  
  # Prepare parasitoïde data
  heatmap_para_melted <- melt(heatmap_densite_para, id.vars = "Row.Labels")
  colnames(heatmap_para_melted) <- c("Y", "X", "Density")
  heatmap_para_melted$X <- as.numeric(as.character(heatmap_para_melted$X))
  
  # Prepare limit existence 
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
  
  
  # Create the base heatmap plot
  axeY <- paste0(" mortality rate of ",variable_legende)
  
  #Preparer les données transformé
  heatmap_plot <- heatmap_para_melted
  
  for (i in 1:length(heatmap_plot$Density)) {
    if (heatmap_melted$Density[i] < 1) {
      heatmap_plot$Density[i] = "All dead"
    } else {
      if (heatmap_para_melted$Density[i] < 0.01){
        if (isocline_melted$Density[i] < 0.01 ) {
          heatmap_plot$Density[i] = "Stable without parasitoid"
        } else {
          heatmap_plot$Density[i] = "Unstable without parasitoid"
        }
      } else {
        if (isocline_melted$Density[i] < 0.01 ) {
          heatmap_plot$Density[i] = "Stable with parasitoid"
        } else {
          heatmap_plot$Density[i] = "Unstable with parasitoid"
        }
      }
      
    }
  }
  
  heatmap_plot$Density <- factor(heatmap_plot$Density, levels = c("All dead", "Stable without parasitoid", "Unstable without parasitoid", "Stable with parasitoid", "Unstable with parasitoid"))
  pastel_colors <- brewer.pal(n = 5, name = "Pastel1")
  # Ta HeatMap avec couleurs définies manuellement
  p <- ggplot(heatmap_plot) +
    geom_tile(aes(x = X, y = Y, fill = Density)) +
    scale_fill_manual(
      name = "System state",
      #values = c(
      #  "All dead"= pal_material("grey")(10)[7],
      #  "Stable without parasitoid" = pal_material("light-blue")(10)[5],
      #  "Unstable without parasitoid" = pal_material("light-blue")(10)[7],
      #  "Stable with parasitoid" = pal_material("light-green")(10)[5],
      #  "Unstable with parasitoid" = pal_material("light-green")(10)[7]
      #)
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
  p = p + geom_line(data = isocline_existance, aes(x = X, y = Y), 
                    color = "#666666", linetype = "dashed", linewidth = 1)
  
  p <- p +
    theme(
      panel.background = element_rect(fill = "#FFF6E5" , colour = NA),
      plot.background  = element_rect(fill = "#FFF6E5" , colour = NA)
    )
  
  return(p)
}


# Base function for the creation HeatMap states
create_plot_Relatif <- function(pop_obs ,model_name, variable_name, ligne_name, variable_legende, isDD, racine, isLimite, global_min, global_max, Xtaille, x_grille, y_grille, title_names) {
  # Load the data
  # Densité Oeufs
  #file_path <- paste0(racine, model_name, '_heatmap_Densité_',isDD,'_',variable_name[2],'_', variable_name[1], '.csv')
  
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


# Function to find global min and max for heatmap data
find_global_min_max <- function(pop_obs ,racine, model_names, variable_names1, variable_names2, isDD) {
  min_val <- Inf
  max_val <- -Inf
  for (model_name in model_names) {
    for (variable_name1 in variable_names1) {
      for (variable_name2 in variable_names2) {
        # Densité des ravageur classique
        file_path <- paste0(racine, model_name, pop_obs, isDD,"_", variable_name2,"_", variable_name1, '.csv')
        
        if (file.exists(file_path)) {
          heatmap_densite <- read.csv(file_path)
          heatmap_densite <- heatmap_densite[, -1]  # Remove row labels
          heatmap_densite <- as.matrix(heatmap_densite)
          heatmap_densite <- as.numeric(heatmap_densite)
          
          current_min <- min(heatmap_densite, na.rm = TRUE)
          current_max <- max(heatmap_densite, na.rm = TRUE)
          
          if (current_min < min_val) min_val <- current_min
          if (current_max > max_val) max_val <- current_max
          
          #if (min_val == 0) min_val <- 0.01
        }
      }
    }
  }
  
  return(list(min = min_val, max = max_val))
}


# Function to save a plot
save_plot <- function(plot, model_name, isDD, variable_name, racineSave) {
  plot_file_path <- paste0(racineSave, model_name,"_",isDD,"_q_",variable_name, "_heatmap_plot.png")
  ggsave(plot_file_path, plot = plot, width = 10, height = 8)
  cat("Plot saved as", plot_file_path, "\n")
}


# Function to create a model-specific plot
create_model_specific_plot <- function(model_name, Empfile, isLimite) {
  # Check the operating system and set the root directory accordingly
  if (Sys.info()["sysname"] == "Windows") {
    racine <-  paste0("P:/Article1/tables/",Empfile,"/")
  } else {
    racine <- paste0("/Users/Quentin/Desktop/Article1/tables/",Empfile,"/")
  }
  
  variable_names = c("deltaA","deltaU","deltaL","deltaM")
  ligne_names = c("DeltaA1","DeltaU1","DeltaL2","DeltaM2")
  variable_legende = c("ravageurs adultes", "oeufs de ravageur", "larves de ravageur", "nymphe de ravageur")
  
  plots <- list()
  
  for (i in 1:4) {
    variable_name = variable_names[i]
    ligne_name = ligne_names[i]
    
    if (model_name != "sansPara") {
      # Generate plots for both DD and ssDD
      plots[[paste0(variable_name, "_DD")]] <- create_plot(model_name, variable_name, ligne_name, variable_legende[i], "DD", racine, isLimite)
      plots[[paste0(variable_name, "_ssDD")]] <- create_plot(model_name, variable_name, ligne_name, variable_legende[i], "ssDD", racine, isLimite)
    } else {
      # Generate plot for DD only
      plots[[paste0(variable_name, "_DD")]] <- create_plot(model_name, variable_name, ligne_name, variable_legende[i], "DD", racine, isLimite)
      plots[[paste0(variable_name, "_ssDD")]] <- NULL  # Empty plot for sansPara
    }
  }
  
  # Arrange plots in a grid
  grid_plot <- grid.arrange(
    arrangeGrob(grobs = plots[grepl("_DD", names(plots))], ncol = 1),
    arrangeGrob(grobs = plots[grepl("_ssDD", names(plots))], ncol = 1),
    ncol = 2
  )
  
  # Save the combined plot
  if (Sys.info()["sysname"] == "Windows") {
    racineSave <-  paste0("P:/Article1/tables/graphs/graphR/",Empfile,"/")
  } else {
    racineSave <- paste0("/Users/Quentin/Desktop/Article1/graphs/graphR/",Empfile,"/")
  }
  plot_file_path <- paste0(racineSave, model_name, "_specific_plot.png")
  ggsave(plot_file_path, plot = grid_plot, width = 20, height = 24)
  cat("Model-specific plot saved as", plot_file_path, "\n")
}


# Function to create comparison tables
create_comparison_tables <- function(Empfile, isLimite, type = "Densité",variable_names = c("deltaU","deltaL","deltaM","deltaA"),variable2_names="q",ligne_names = c("DeltaA1","DeltaU1","DeltaL2","DeltaM2"),
                                     variable_legende = c("egg", "larva", "pupa", "adult"), title_names = c("Without parasitoid","Egg parasitoid","Larva parasitoid","Pupa parasitoid")) {
  # Check the operating system and set the root directory accordingly
  if (Sys.info()["sysname"] == "Windows") {
    racine <-  paste0("P:/Article1/tables/",Empfile,"/")
  } else {
    racine <- paste0("/Users/Quentin/Desktop/Article1/tables/",Empfile,"/")
  }
  #pops_obs <- c('_heatmap_Densité_Oeuf_','_heatmap_Densité_','_heatmap_Densité_Nymphe_','_heatmap_Densité_Adulte_','_heatmap_Densité_ParaJuv_','_heatmap_Densité_Para_')
  pops_obs <- c('_heatmap_Densité_')
  model_names = c("sansPara", "ParaOeufs", "ParaLarve", "ParaNymphe")
  
  
  for (pop_obs in pops_obs) {
    # Find global min and max for DD and ssDD
    global_min_max_DD <- find_global_min_max(pop_obs, racine, model_names, variable_names[1:3], variable2_names,"DD")
    global_min_max_ssDD <- find_global_min_max(pop_obs, racine, model_names[-1], variable_names[1:3], variable2_names, "ssDD")
    
    
    
    plots_DD <- list()
    plots_ssDD <- list()
    
    o = 0
    count_col = 1
    for (i in 1:length(variable_names)) {
      variable_name = variable_names[i]
      ligne_name = ligne_names[i]
      o = o+1
      model_name = model_names[1]
      if (type == "Densité" && variable2_names == "q"){
        # Generate plot for DD only
        plots_DD[[paste0(model_name, "_", variable_name)]] <- create_plot(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                          "DD", racine, isLimite, global_min_max_DD$min, global_min_max_DD$max, 70, o, i, title_names[count_col])
        #plots_ssDD[[paste0(model_name, "_", variable_name)]] <- NULL  # Empty plot for sansPara
      } else if (type == "Type"&& variable2_names == "q") {
        plots_DD[[paste0(model_name, "_", variable_name)]] <- create_plot_States(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                                 "DD", racine, isLimite, global_min_max_DD$min, global_min_max_DD$max, 70, o, i, title_names[count_col])
      }
      
      for (model_name in model_names[2:length(model_names)]) {
        o = o+1
        count_col = count_col + 1
        
        if (type == "Densité"){
          # Generate plots for DD
          plots_DD[[paste0(model_name, "_", variable_name)]] <- create_plot(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                            "DD", racine, isLimite, global_min_max_DD$min, global_min_max_DD$max, 70, o, i, title_names[count_col])
          #plots_ssDD[[paste0(model_name, "_", variable_name)]] <- create_plot(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
          #                                                                    "ssDD", racine, isLimite, global_min_max_DD$min, global_min_max_DD$max, 70, o, i)
        } else if (type == "Type") {
          plots_DD[[paste0(model_name, "_", variable_name)]] <- create_plot_States(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                                   "DD", racine, isLimite, global_min_max_DD$min, global_min_max_DD$max, 70, o, i, title_names[count_col])
        }else if (type == "Relatif") {
          plots_DD[[paste0(model_name, "_", variable_name)]] <- create_plot_Relatif(pop_obs, model_name, c(variable_name,variable2_names), ligne_name, variable_legende[i], 
                                                                                    "DD", racine, isLimite, global_min_max_DD$min, global_min_max_DD$max, 70, o, i, title_names[count_col])
        }
      }
    }
    
    # Save the combined plots
    if (Sys.info()["sysname"] == "Windows") {
      racineSave <-  paste0("P:/Article1/tables/graphs/graphR/",Empfile,"/")
    } else {
      racineSave <- paste0("/Users/Quentin/Desktop/Article1/graphs/graphR/",Empfile,"/")
    }
    plot_file_path_DD <- paste0(racineSave, pop_obs, "_", type, "_", "comparison_table_DD.png")
    #plot_file_path_ssDD <- paste0(racineSave,  pop_obs, "_", type, "_", "comparison_table_ssDD.png")
    
    # Arrange plots in tables
    if (type == "Relatif" ||variable2_names!="q"){
      grid_plot_DD <- grid.arrange(grobs = plots_DD, ncol = 3)
      ggsave(plot_file_path_DD, plot = grid_plot_DD, width = 20, height = 20, dpi = 300)
      #ggsave(plot_file_path_ssDD, plot = grid_plot_ssDD, width = 30, height = 20)
    }else if (variable_names[1] !="alphaL" && length(variable_names)==1){
      grid_plot_DD <- grid.arrange(grobs = plots_DD, ncol = 3)
      ggsave(plot_file_path_DD, plot = grid_plot_DD, width = 20, height = 20, dpi = 300)
    }else{
      #plots_DD_clean <- lapply(plots_DD, function(p) p + theme(legend.position = "none"))
      #grid_plot_DD <- grid.arrange(grobs = plots_DD_clean, ncol = 4)
      ##grid_plot_ssDD <- grid.arrange(grobs = plots_ssDD, ncol = 4)
      #ggsave(plot_file_path_DD, plot = grid_plot_DD, width = 30, height = 20, dpi = 300)
      ##ggsave(plot_file_path_ssDD, plot = grid_plot_ssDD, width = 30, height = 20)
      
      
      
      # 1. Récupérer la légende depuis un des plots d'origine
      legend_plot <- plots_DD[[10]] + 
        theme(
          legend.position = "right",
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 30, face = "bold"),
          legend.key.size = unit(1.2, "cm"),
          legend.spacing.y = unit(0.5, "cm"),
          panel.background  = element_rect(fill = hcl(39, 100, 95), colour = NA),
          plot.background   = element_rect(fill = hcl(39, 100, 95), colour = NA)
        )
      legend <- cowplot::get_legend(legend_plot )
      
      # 2. Convertir en un ggdraw pour l'utiliser dans patchwork
      legend_gg <- cowplot::ggdraw(legend)
      
      
      # Création d'un objet de même dimension que la kégende mais vide
      empty_legend <- plots_DD[[10]] + 
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
      
      # 3. Nettoyer tous les plots
      plots_DD_clean <- lapply(plots_DD, function(p) {
        p + 
          theme(
            legend.position = "none",
            plot.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()
          )
      })
      
      # 4. Supprimer axes internes
      #for (i in 1:16) {
      #  if (i <= 12) {
      #    plots_DD_clean[[i]] <- plots_DD_clean[[i]] +
      #      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      #  }
      #  if (i %% 4 != 1) {
      #    plots_DD_clean[[i]] <- plots_DD_clean[[i]] +
      #      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      #  }
      #}
      
      # 5. Créer titres colonnes et lignes
      col_titles <- c("Without parasitoid", "Egg parasitoid", "Larva parasitoid", "Pupa parasitoid")
      row_labels <- c("Egg mortality rate", "Larva mortality rate", "Pupa mortality rate", "Adult mortality rate")
      
      # 6. Créer chaque rangée de graphiques avec titre de ligne
      plot_rows <- list()
      for (i in 0:3) {
        label_grob <- ggplot() +
          theme_void() +
          annotation_custom(
            grob = textGrob(row_labels[i + 1], rot = 90, gp = gpar(fontsize = 30)),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
          )
        row_plots <- plots_DD_clean[(1 + i*4):(4 + i*4)]
        
        # La légende uniquement pour la 1ère ligne
        legend_gg <- if (i == 0) legend_gg else empty_legend_gg
        
        row <- label_grob + row_plots[[1]] + row_plots[[2]] + row_plots[[3]] + row_plots[[4]] + legend_gg
        plot_rows[[i + 1]] <- row + plot_layout(ncol = 6, widths = c(0.1, rep(1, 4), 1.1)) # 7 colonne, 4 de 1 de large (les schéma), une de 0.1 et une de 0.2 (la légende)
      }
      
      # 7. Ligne du haut avec titres de colonnes
      col_title_plots <- lapply(col_titles, function(title) {
        ggplot() +
          theme_void() +
          annotation_custom(textGrob(title, gp = gpar(fontsize = 30, fontface = "bold")))
      })
      empty_plot <- ggplot() + theme_void()
      empty_legend_spot <- ggplot() + theme_void()  # case vide pour la légende
      top_row <- empty_plot + col_title_plots[[1]] + col_title_plots[[2]] + col_title_plots[[3]] + col_title_plots[[4]] + empty_legend_spot 
      top_row <- top_row + plot_layout(ncol = 6, widths = c(0.1, rep(1, 4), 1.1))
      
      # 8. Axe X commun
      bottom_axis <- ggplot() +
        theme_void() +
        annotation_custom(
          grob = textGrob("Fecondity rate", gp = gpar(fontsize = 30))
        )
      
      # 9. Fusion finale
      final_plot <- top_row /
        plot_spacer() /
        plot_rows[[1]] /
        plot_rows[[2]] /
        plot_rows[[3]] /
        plot_rows[[4]] /
        bottom_axis
      
      final_plot <- final_plot +
        plot_layout(heights = c(0.1, 0.01, rep(1, 4), 0.1))
      
      # 10. Affichage
      print(final_plot)
      
      # 11. Sauvegarde
      
      ggsave(plot_file_path_DD, plot = final_plot, width = 25, height = 23, dpi = 300,   bg ="#FFF6E5")
    }
    
    
    cat("Comparison tables saved as", plot_file_path_DD, "and", "plot_file_path_ssDD", "\n")
    
    return()
  }
  
}











#######################################
########### Figure à tracer ###########
#######################################


## Final pour article 

### HeatMap densité Final (dev des para 16)
a=create_comparison_tables("HeatMapDD_paradev_16_final/fixeSimple/0.0001",FALSE)

a
### HeatMap densité Relative Final (dev des para 16)
create_comparison_tables("HeatMapDD_paradev_16_final/fixeSimple/0.0001",FALSE, "Relatif")

### HeatMap densité Final (dev des para 16) version type d'équilibre
create_comparison_tables("HeatMapDD_paradev_16_final/fixeSimple/0.0001",FALSE, "Type")

### HeatMap densité Final (dev des para 16) version type d'équilibre
create_comparison_tables("HeatMapDD_paradev_DevMorta/fixeSimple/0.0001",FALSE, variable_names = c("deltaU","deltaL","deltaM","deltaA"),variable2_names="tauJ",
                         ligne_names = c("DeltaA1","DeltaU1","DeltaL2","DeltaM2"), variable_legende = c("oeufs", "larves", "nymphes","adultes"))











## La Grande modification 


#### taux d'attaque 

###### a = 0
create_comparison_tables("GrandeVerification/Verif_heatmap_a0/fixeSimple/0.0001",FALSE, "Type")
###### a = 0.03
create_comparison_tables("GrandeVerification/Verif_heatmap_a003/fixeSimple/0.0001",FALSE, "Type")
###### a = 0.001
create_comparison_tables("GrandeVerification/Verif_heatmap_a0001/fixeSimple/0.0001",FALSE, "Type")



#### Durée de vie des stades juv de 4
create_comparison_tables("GrandeVerification/Verif_heatmap_tpsdev4/fixeSimple/0.0001",FALSE, "Type")



#### mortalité de base des stades juv de 0.1
create_comparison_tables("GrandeVerification/Verif_heatmap_delta01/fixeSimple/0.0001",FALSE, "Type")



#### Variation de la compétition 

###### alphaL = 0.00001
create_comparison_tables("GrandeVerification/Verif_heatmap_alphaL000001/fixeSimple/0.0001",FALSE, "Type")
###### alphaL = 0.0000001
create_comparison_tables("GrandeVerification/Verif_heatmap_alphaL00000001/fixeSimple/0.0001",FALSE, "Type")



#### Exposant de la competition #pas calculer jusu'au bout, regarder les LOG

###### exposant = 1.1
create_comparison_tables("GrandeVerification/Verif_heatmap_DD11/fixeSimple/0.0001",FALSE, "Type")
###### exposant = 1.5
create_comparison_tables("GrandeVerification/Verif_heatmap_DD15/fixeSimple/0.0001",FALSE, "Type")
###### exposant = 1.8
create_comparison_tables("GrandeVerification/Verif_heatmap_DD18/fixeSimple/0.0001",FALSE, "Type")



#### Variation du temps de dev des ravageurs 

###### tauL = 1
create_comparison_tables("GrandeVerification/Verif_heatmap_tauJ1/fixeSimple/0.0001",FALSE, "Type")
###### tauL = 14
create_comparison_tables("GrandeVerification/Verif_heatmap_tauJ14/fixeSimple/0.0001",FALSE, "Type")
###### tauL = 28
create_comparison_tables("GrandeVerification/Verif_heatmap_tauJ28/fixeSimple/0.0001",FALSE, "Type")


#### RF de type II : variation du temps de manipulation

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



