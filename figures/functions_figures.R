library(ggplot2)
library(ggforce)
# Load necessary libraries
library(spectralBridges)
library(ggplot2)
library(ClusterR)
library(deldir)  # for computing Voronoi cells
library(ggforce) # for drawing the Voronoi cells

set.seed(0)
nk <- 30
mu2 <- c(4, 4)
X1 <- matrix(rnorm(nk * 2), nk, 2)
X2 <- matrix(rnorm(nk * 2), nk, 2) + matrix(mu2, nk, 2, byrow = TRUE)
clusters <- rep(c(1, 2), each = nk)
center <- mu2 / 2

plot_balls <- function(a, b, c, d, xlim = c(-5, 7), ylim = c(-5, 7)) {
  # Calculer les coefficients de la droite orthogonale au vecteur (c, d)
  slope <- -c / d
  
  # Définir la fonction de la droite
  line_function <- function(x) {
    return(slope * (x - a) + b)
  }
  
  # Créer une séquence de points x pour tracer la droite
  x_values <- seq(xlim[1], xlim[2], length.out = 100)
  y_values <- line_function(x_values)
  
  # Créer des data frames pour ggplot
  df_line <- data.frame(x = x_values, y = y_values)
  df_point <- data.frame(x = a, y = b)
  df_clusters <- data.frame(x = rbind(X1, X2)[,1], y = rbind(X1, X2)[,2], cluster = factor(clusters))
  df_vector <- data.frame(x = c(a, a + c), y = c(b, b + d))
  df_centers <- data.frame(x = c(0, mu2[1]), y = c(0, mu2[2]), label = c("Centroid k", "Centroid l"))
  df_segments <- rbind(
    data.frame(x = X1[, 1], y = X1[, 2], xend = 0, yend = 0),
    data.frame(x = X2[, 1], y = X2[, 2], xend = mu2[1], yend = mu2[2])
  )
  
  # Tracer avec ggplot2
  p <- ggplot() +
    geom_segment(data = df_segments, aes(x = x, y = y, xend = xend, yend = yend), color = 'gray', size = 0.8) +
    geom_line(data = df_line, aes(x = x, y = y), color = 'black') +
    geom_point(data = df_clusters, aes(x = x, y = y, color = cluster)) +
    geom_point(data = df_centers, aes(x = x, y = y), color = 'black', size = 3) +
    geom_ellipse(aes(x0 = 0, y0 = 0, a = 2.5, b = 2.5, angle = 0), color = 'red', linetype = 2, fill = 'red', alpha = 0.2) +
    geom_ellipse(aes(x0 = mu2[1], y0 = mu2[2], a = 2.5, b = 2.5, angle = 0), color = 'green', linetype = 2, fill = 'green', alpha = 0.2) +
    xlim(xlim) +
    ylim(ylim) +
    theme_minimal() +
    theme(legend.position = "topright") +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )
  
  print(p)
}


# Function to compute the projection of point x onto segment ab
projection <- function(x, a, b) {
  ab <- b - a
  t <- sum((x - a) * ab) / sum(ab * ab)
  t <- max(0, min(1, t))
  a + t * ab
}

plot_bridge <- function(a, b, c, d, xlim = c(-5, 7), ylim = c(-5, 7)) {
  # Calculate the slope of the line orthogonal to vector (c, d)
  slope <- -c / d
  
  # Define the function of the line
  line_function <- function(x) {
    return(slope * (x - a) + b)
  }
  
  # Create a sequence of x points to plot the line
  x_values <- seq(xlim[1], xlim[2], length.out = 100)
  y_values <- line_function(x_values)
  
  # Create data frames for ggplot
  df_line <- data.frame(x = x_values, y = y_values)
  df_clusters <- data.frame(x = rbind(X1, X2)[, 1], y = rbind(X1, X2)[, 2], cluster = factor(clusters))
  df_centers <- data.frame(x = c(0, mu2[1]), y = c(0, mu2[2]), label = c("Centroid k", "Centroid l"))
  
  # Calculer la longueur et le centre du rectangle
  rect_length <- sqrt(sum((mu2 - c(0, 0))^2))
  rect_center <- (mu2 + c(0, 0)) / 2
  rect_angle <- atan2(mu2[2], mu2[1])
  
  # Calculer les coins du rectangle
  rect_width <- 5
  half_length <- rect_length / 2
  half_width <- rect_width / 2
  
  x1 <- rect_center[1] - half_length * cos(rect_angle) - half_width * sin(rect_angle)
  y1 <- rect_center[2] - half_length * sin(rect_angle) + half_width * cos(rect_angle)
  x2 <- rect_center[1] + half_length * cos(rect_angle) - half_width * sin(rect_angle)
  y2 <- rect_center[2] + half_length * sin(rect_angle) + half_width * cos(rect_angle)
  x3 <- rect_center[1] + half_length * cos(rect_angle) + half_width * sin(rect_angle)
  y3 <- rect_center[2] + half_length * sin(rect_angle) - half_width * cos(rect_angle)
  x4 <- rect_center[1] - half_length * cos(rect_angle) + half_width * sin(rect_angle)
  y4 <- rect_center[2] - half_length * sin(rect_angle) - half_width * cos(rect_angle)
  
  df_rect <- data.frame(
    x = c(x1, x2, x3, x4),
    y = c(y1, y2, y3, y4),
    id = rep(1, 4)
  )
  
  
  
  
  # Calculate projections and segments
  a <- c(0, 0)
  b <- mu2
  projections <- apply(rbind(X1, X2), 1, projection, a, b)
  projections <- t(projections)
  
  segments <- data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
  
  for (i in 1:nrow(projections)) {
    point <- rbind(X1, X2)[i, ]
    proj <- projections[i, ]
    centroid <- if (clusters[i] == 1) a else b
    if (sum((point - proj)^2) < sum((point - centroid)^2)) {
      segments <- rbind(segments, data.frame(x = point[1], y = point[2], xend = proj[1], yend = proj[2]))
    } else {
      segments <- rbind(segments, data.frame(x = point[1], y = point[2], xend = centroid[1], yend = centroid[2]))
    }
  }
  
  # Plot with ggplot2
  p <- ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = mu2[1], yend = mu2[2]), color = 'blue') + # Segment joining centroids
    geom_segment(data = segments, aes(x = x, y = y, xend = xend, yend = yend), color = 'gray', size = 0.8) +
    geom_line(data = df_line, aes(x = x, y = y), color = 'black') +
    geom_point(data = df_clusters, aes(x = x, y = y, color = cluster)) +
    geom_point(data = df_centers, aes(x = x, y = y), color = 'black', size = 3) +
    geom_ellipse(aes(x0 = 0, y0 = 0, a = 2.5, b = 2.5, angle = 0), color = 'red', linetype = 2, fill = 'red', alpha = 0.2) +
    geom_ellipse(aes(x0 = mu2[1], y0 = mu2[2], a = 2.5, b = 2.5, angle = 0), color = 'green', linetype = 2, fill = 'green', alpha = 0.2) +
    geom_polygon(data = df_rect, aes(x = x, y = y, group = id), fill = "blue", alpha = 0.2) +
    xlim(xlim) +
    ylim(ylim) +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )
  
  print(p)
}



#pdf("balls.pdf")
plot_balls(a = center[1], b = center[2], c = mu2[1], d = mu2[2])
#dev.off()

# Exemple d'utilisation
# pdf("bidge.pdf")
plot_bridge(a = center[1], b = center[2], c = mu2[1], d = mu2[2])
# dev.off()




# Load the iris dataset
data(iris)

# Perform PCA on the iris dataset
pca <- prcomp(iris[, 1:4], scale. = TRUE)
pca_data <- as.data.frame(pca$x)
X<-pca_data[, 1:2]

# Use KMeans_rcpp to cluster the PCA-transformed data into 20 clusters
set.seed(123)  # for reproducibility
n_nodes=20
n<-nrow(X)
kmeans_result  = KMeans_rcpp(X, clusters = n_nodes, 
                             num_init = 3, max_iters =30, 
                             initializer = 'kmeans++')

kmeans_centers <- as.matrix(kmeans_result$centroids)
kmeans_labels <-  as.matrix(kmeans_result$clusters)
kmeans_size <-    kmeans_result$obs_per_cluster
kmeans_Iw<-kmeans_result$WCSS_per_cluster

# 2. Affinity computation 
###################################

affinity <- matrix(0, n_nodes, n_nodes)
X.centered <-
  as.matrix(do.call(rbind, lapply(1:n, function(i) {
    X[i, ] - kmeans_centers[kmeans_labels[i], ]
  })))

affinity <- matrix(0, n_nodes, n_nodes)
X.centered <-
  as.matrix(do.call(rbind, lapply(1:n, function(i) {
    X[i, ] - kmeans_centers[kmeans_labels[i], ]
  })))

affinity <- outer(1:n_nodes, 1:n_nodes,
                  Vectorize(function(k, l) {
                    if (k==l) 0 else {
                      distkl2<-sum(((kmeans_centers[l,] - kmeans_centers[k, ]))^2);        
                      alphai<-c(pmax(0, (kmeans_centers[l,] - kmeans_centers[k, ]) %*% 
                                       t(X.centered[(kmeans_labels ==k), ]))/distkl2,
                                (pmax(0, (kmeans_centers[k,] - kmeans_centers[l, ]) %*% 
                                        t(X.centered[(kmeans_labels ==l), ]))/distkl2));
                      sqrt(sum(alphai^2)/(kmeans_size[k]+kmeans_size[l]))
                      #compute_entropy(ti)
                    }
                  }))




# Get the cluster centroids
centroids <- as.data.frame(kmeans_result$centroids)
colnames(centroids) <- c("PC1", "PC2")

# Compute Voronoi cells
voronoi <- deldir(centroids$PC1, centroids$PC2)

# Extract Voronoi edges
voronoi_edges <- voronoi$dirsgs

# Convert edges to a data frame for ggplot2
voronoi_edges_df <- as.data.frame(voronoi_edges)

# Plot the first two PCA components with Voronoi cells
p <- ggplot() +
  geom_point(data = X, aes(x = PC1, y = PC2, color = as.factor(kmeans_result$clusters),alpha=0.4)) +
  geom_polygon(data = voronoi_edges_df, aes(x = x, y = y, group = id), color = "black", alpha = 0.2,fill="white") +
  geom_point(data = centroids, aes(x = PC1, y = PC2), color = "red", size = 5) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
print(p)

