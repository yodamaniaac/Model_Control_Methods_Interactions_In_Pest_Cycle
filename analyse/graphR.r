# Load necessary libraries
library(ggplot2)
library(reshape2)
library(viridis)


# Check the operating system and set the root directory accordingly
if (Sys.info()["sysname"] == "Windows") {
  racine <- "P:/Article1/MeSU/REZ/Model0/"
} else {
  racine <- "/Users/Quentin/Desktop/Article1/MeSU/REZ/Model0/"
}

# Load the data
file_path <- paste0(racine,'sansPara_heatmap_Densité_DD_q_deltaA.csv')
heatmap_densite <- read.csv(file_path)
colnames(heatmap_densite) <- sub("^X", "", colnames(heatmap_densite))
head(heatmap_densite)

file_path <- paste0(racine,'sansPara_heatmap_Variance_DD_q_deltaA.csv')
heatmap_Variance <- read.csv(file_path)
colnames(heatmap_Variance) <- sub("^X", "", colnames(heatmap_Variance))
head(heatmap_Variance)

file_path <- paste0(racine,'LimiteStabDeltaA1.csv')
heatmap_Limite <- read.csv(file_path, header = TRUE)
head(heatmap_Limite)
colnames(heatmap_Limite) <- sub("^X", "", colnames(heatmap_Limite))
heatmap_Limite <- rbind(colnames(heatmap_Limite), heatmap_Limite)
colnames(heatmap_Limite) <- c("X", "Y")
head(heatmap_Limite)
heatmap_Limite$X <- as.numeric(as.character(heatmap_Limite$X))
heatmap_Limite$Y <- as.numeric(as.character(heatmap_Limite$Y))
heatmap_Limite <- heatmap_Limite[heatmap_Limite$X <= 70, ]
head(heatmap_Limite)
max(heatmap_Limite$X)
length(heatmap_Limite$X)



# Prepare heatmap data
heatmap_melted <- melt(heatmap_densite, id.vars = "Row.Labels")
colnames(heatmap_melted) <- c("Y", "X", "Density")
heatmap_melted$X <- as.numeric(as.character(heatmap_melted$X))
head(heatmap_melted)

# Prepare isocline data
isocline_melted <- melt(heatmap_Variance, id.vars = "Row.Labels")
colnames(isocline_melted) <- c("Y", "X", "Density")
isocline_melted$X <- as.numeric(as.character(isocline_melted$X))
head(isocline_melted)

# Create the base heatmap plot
p <- ggplot(heatmap_melted) +
  geom_tile(aes(x = X, y = Y, fill = Density)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  scale_x_continuous(expand = c(0, 0)) +  # Supprime l'espace autour des cases
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 70) +
  labs(title = "Heatmap with Isoclines and Stability Curve",
       x = "X-axis Label",
       y = "Y-axis Label")

print(p)

# Add isoclines
q <- p + geom_contour(data = isocline_melted, aes(x = X, y = Y, z = Density),
                      color = "red", alpha =1, breaks = c(1,100,200,300), linewidth = 1)

print(q)

# Add stability limit curve
heatmap_Limite <- heatmap_Limite[order(heatmap_Limite$Y), ]
r <- q + geom_point(data = heatmap_Limite, aes( y = Y,x = X), color = "grey", size = 0.01, alpha = 10) 
print(r)

# Display the plot
print(r)

