# Load necessary libraries
library(ggplot2)
library(reshape2)
library(viridis)

# Define the function
create_and_save_plot <- function(model_name) {
  # Check the operating system and set the root directory accordingly
  if (Sys.info()["sysname"] == "Windows") {
    racine <- "P:/Article1/MeSU/REZ/Model0/"
  } else {
    racine <- "/Users/Quentin/Desktop/Article1/MeSU/REZ/Model0/"
  }
  
  variable_names = c("deltaA","deltaU","deltaL","deltaM")
  ligne_names = c("DeltaA1","DeltaU1","DeltaL2","DeltaM2")
  
  for (i in 1:4) {
  
    #Prendre le nom des variable
    variable_name = variable_names[i]
    ligne_name = ligne_names[i]
    
    # Load the data
    file_path <- paste0(racine, model_name, '_heatmap_Densité_DD_q_', variable_name, '.csv')
    heatmap_densite <- read.csv(file_path)
    colnames(heatmap_densite) <- sub("^X", "", colnames(heatmap_densite))
    
    file_path <- paste0(racine, model_name, '_heatmap_Variance_DD_q_', variable_name, '.csv')
    heatmap_Variance <- read.csv(file_path)
    colnames(heatmap_Variance) <- sub("^X", "", colnames(heatmap_Variance))
    
    file_path <- paste0(racine, 'LimiteStab',ligne_name,'.csv')
    heatmap_Limite <- read.csv(file_path, header = TRUE)
    colnames(heatmap_Limite) <- sub("^X", "", colnames(heatmap_Limite))
    heatmap_Limite <- rbind(colnames(heatmap_Limite), heatmap_Limite)
    colnames(heatmap_Limite) <- c("X", "Y")
    heatmap_Limite$X <- as.numeric(as.character(heatmap_Limite$X))
    heatmap_Limite$Y <- as.numeric(as.character(heatmap_Limite$Y))
    heatmap_Limite <- heatmap_Limite[heatmap_Limite$X <= 70, ]
    
    # Prepare heatmap data
    heatmap_melted <- melt(heatmap_densite, id.vars = "Row.Labels")
    colnames(heatmap_melted) <- c("Y", "X", "Density")
    heatmap_melted$X <- as.numeric(as.character(heatmap_melted$X))
    
    # Prepare isocline data
    isocline_melted <- melt(heatmap_Variance, id.vars = "Row.Labels")
    colnames(isocline_melted) <- c("Y", "X", "Density")
    isocline_melted$X <- as.numeric(as.character(isocline_melted$X))
    
    # Create the base heatmap plot
    p <- ggplot(heatmap_melted) +
      geom_tile(aes(x = X, y = Y, fill = Density)) +
      scale_fill_viridis_c() +
      theme_minimal() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed(ratio = 70) +
      labs(title = "Heatmap with Isoclines and Stability Curve",
           x = "X-axis Label",
           y = "Y-axis Label")
    
    # Add isoclines
    q <- p + geom_contour(data = isocline_melted, aes(x = X, y = Y, z = Density),
                          color ="yellow", alpha = 1, breaks = c(1), linewidth = 1) +
             geom_contour(data = isocline_melted, aes(x = X, y = Y, z = Density),
                          color ="orange", alpha = 1, breaks = c(100), linewidth = 1) +
             geom_contour(data = isocline_melted, aes(x = X, y = Y, z = Density),
                          color ="red", alpha = 1, breaks = c(200), linewidth = 1)
    
    # Add stability limit curve
    heatmap_Limite <- heatmap_Limite[order(heatmap_Limite$Y), ]
    r <- q + geom_point(data = heatmap_Limite, aes(y = Y, x = X), color = "white", size = 0.01, alpha = 10)
    
    # Save the plot
    plot_file_path <- paste0(racine, model_name,"_q_",variable_name, "_heatmap_plot.png")
    ggsave(plot_file_path, plot = r, width = 10, height = 8)
    print(r)
    # Print a message indicating the plot has been saved
    cat("Plot saved as", plot_file_path, "\n")
  }
  
}

# Example usage
create_and_save_plot('sansPara')

