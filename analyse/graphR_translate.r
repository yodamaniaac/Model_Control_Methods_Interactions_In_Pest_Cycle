# Load necessary libraries
library(ggplot2)
library(reshape2)


# Check the operating system and set the root directory accordingly
if (Sys.info()["sysname"] == "Windows") {
  racine <- "P:/Article1/MeSU/REZ/Model0/"
} else {
  racine <- "/Users/Quentin/Desktop/Article1/MeSU/REZ/Model0/"
}

# Load the data
file_path <- paste0(racine,'sansPara_heatmap_Densité_DD_q_deltaA.csv')
heatmap_densite <- read.csv(file_path)
row_names <- heatmap_densite$Row.Labels
heatmap_densite$Row.Labels <- NULL
colnames(heatmap_densite) <- sub("^X", "", colnames(heatmap_densite))
col_names <- colnames(heatmap_densite)

# Convertir le dataframe en matrice
heatmap_densite_matrix <- as.matrix(heatmap_densite)
colnames(heatmap_densite_matrix) <- NULL

# Permuter les dimensions
heatmap_densite_permuted <- aperm(heatmap_densite_matrix, perm = c(2, 1))

# Repasser en dataframe
heatmap_densite <- as.data.frame(heatmap_densite_permuted)
rownames(heatmap_densite) <- row_names
colnames(heatmap_densite) <- col_names
heatmap_densite$Row.Labels <- row_names
head(heatmap_densite)




# Load the data
file_path <- paste0(racine,'sansPara_heatmap_Variance_DD_q_deltaA.csv')
heatmap_Variance <- read.csv(file_path)
row_names <- heatmap_Variance$Row.Labels
heatmap_Variance$Row.Labels <- NULL
colnames(heatmap_Variance) <- sub("^X", "", colnames(heatmap_Variance))
col_names <- colnames(heatmap_Variance)

# Convertir le dataframe en matrice
heatmap_Variance_matrix <- as.matrix(heatmap_Variance)
colnames(heatmap_Variance_matrix) <- NULL
head(heatmap_Variance_matrix)

# Permuter les dimensions
heatmap_Variance_permuted <- aperm(heatmap_Variance_matrix, perm = c(2, 1))
head(heatmap_Variance_permuted)

# Repasser en dataframe
heatmap_Variance <- as.data.frame(heatmap_Variance_permuted)
rownames(heatmap_Variance) <- row_names
colnames(heatmap_Variance) <- col_names
heatmap_Variance$Row.Labels <- row_names
head(heatmap_densite)


file_path <- paste0(racine,'LimiteStabDeltaA2.csv')
heatmap_Limite <- read.csv(file_path)
colnames(heatmap_Limite) <- sub("^X", "", colnames(heatmap_Limite))
head(heatmap_Limite)



# Prepare heatmap data
heatmap_melted <- melt(heatmap_densite, id.vars = "Row.Labels")
colnames(heatmap_melted) <- c("Y", "X", "Density")
heatmap_melted$X <- as.numeric(as.character(heatmap_melted$X))
head(heatmap_melted)

# Prepare isocline data
isocline_melted <- melt(heatmap_Variance, id.vars = "Row.Labels")
colnames(isocline_melted) <- c("Y", "X", "Variance")
isocline_melted$X <- as.numeric(as.character(isocline_melted$X))
head(isocline_melted)

# Create the base heatmap plot
p <- ggplot(heatmap_melted, aes(x = X, y = Y, fill = Density)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Heatmap with Isoclines and Stability Curve",
       x = "X-axis Label",
       y = "Y-axis Label")

print(p)

# Add isoclines
p <- p + geom_contour(data = isocline_melted, aes(x = X, y = Y, z = Variance),
                      color = "white", alpha = 0.6, linetype = "dashed")

print(p)

# Add stability limit curve
heatmap_Limite$V1 <- as.numeric(as.character(heatmap_Limite$V1))
heatmap_Limite$V2 <- as.numeric(as.character(heatmap_Limite$V2))
p <- p + geom_line(data = heatmap_Limite, aes(x = V1, y = V2), color = "red", size = 1) +
  geom_point(data = heatmap_Limite, aes(x = V1, y = V2), color = "red", size = 1)

print(p)

# Display the plot
print(p)
