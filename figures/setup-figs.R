library(ggplot2)
library(ggforce)
library(ClusterR)
library(deldir)
library(GGally)
library(gridExtra)
library(dplyr)
library(scales)
library(mclust)

set.seed(0)
nk <- 30
mu2 <- c(4, 4)
X1 <- matrix(rnorm(nk * 2), nk, 2)
X2 <- matrix(rnorm(nk * 2), nk, 2) + matrix(mu2, nk, 2, byrow = TRUE)
clusters <- rep(c(1, 2), each = nk)
center <- mu2 / 2

plot_balls <- function(a, b, c, d, xlim = c(-5, 7), ylim = c(-5, 7)) {
  slope <- -c / d
  line_function <- function(x) slope * (x - a) + b
  x_values <- seq(xlim[1], xlim[2], length.out = 100)
  y_values <- line_function(x_values)
  df_line <- data.frame(x = x_values, y = y_values)
  df_clusters <- data.frame(x = rbind(X1, X2)[,1], y = rbind(X1, X2)[,2], cluster = factor(clusters))
  df_centers <- data.frame(x = c(0, mu2[1]), y = c(0, mu2[2]), label = c("Centroid k", "Centroid l"))
  df_segments <- rbind(
    data.frame(x = X1[, 1], y = X1[, 2], xend = 0, yend = 0),
    data.frame(x = X2[, 1], y = X2[, 2], xend = mu2[1], yend = mu2[2])
  )
  p <- ggplot() +
    geom_segment(data = df_segments, aes(x = x, y = y, xend = xend, yend = yend), color = 'gray', linewidth = 0.8) +
    geom_line(data = df_line, aes(x = x, y = y), color = 'black') +
    geom_point(data = df_clusters, aes(x = x, y = y, color = cluster)) +
    geom_point(data = df_centers, aes(x = x, y = y), color = 'black', size = 3) +
    geom_ellipse(aes(x0 = 0, y0 = 0, a = 2.5, b = 2.5, angle = 0), color = 'red', linetype = 2, fill = 'red', alpha = 0.2) +
    geom_ellipse(aes(x0 = mu2[1], y0 = mu2[2], a = 2.5, b = 2.5, angle = 0), color = 'green', linetype = 2, fill = 'green', alpha = 0.2) +
    xlim(xlim) + ylim(ylim) + coord_fixed() +
    theme_minimal() +
    theme(panel.grid = element_blank(), legend.position = "none",
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank())
  p
}

projection <- function(x, a, b) {
  ab <- b - a
  t <- sum((x - a) * ab) / sum(ab * ab)
  t <- max(0, min(1, t))
  a + t * ab
}

plot_bridge <- function(a, b, c, d, xlim = c(-5, 7), ylim = c(-5, 7)) {
  slope <- -c / d
  line_function <- function(x) slope * (x - a) + b
  x_values <- seq(xlim[1], xlim[2], length.out = 100)
  y_values <- line_function(x_values)
  df_line <- data.frame(x = x_values, y = y_values)
  df_clusters <- data.frame(x = rbind(X1, X2)[, 1], y = rbind(X1, X2)[, 2], cluster = factor(clusters))
  df_centers <- data.frame(x = c(0, mu2[1]), y = c(0, mu2[2]), label = c("Centroid k", "Centroid l"))
  rect_length <- sqrt(sum((mu2 - c(0, 0))^2))
  rect_center <- (mu2 + c(0, 0)) / 2
  rect_angle <- atan2(mu2[2], mu2[1])
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
  a0 <- c(0, 0)
  b0 <- mu2
  projections <- t(apply(rbind(X1, X2), 1, projection, a0, b0))
  segments <- data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
  for (i in 1:nrow(projections)) {
    point <- rbind(X1, X2)[i, ]
    proj <- projections[i, ]
    centroid <- if (clusters[i] == 1) a0 else b0
    if (sum((point - proj)^2) < sum((point - centroid)^2)) {
      segments <- rbind(segments, data.frame(x = point[1], y = point[2], xend = proj[1], yend = proj[2]))
    } else {
      segments <- rbind(segments, data.frame(x = point[1], y = point[2], xend = centroid[1], yend = centroid[2]))
    }
  }
  p <- ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = mu2[1], yend = mu2[2]), color = 'blue') +
    geom_segment(data = segments, aes(x = x, y = y, xend = xend, yend = yend), color = 'gray', size = 0.8) +
    geom_line(data = df_line, aes(x = x, y = y), color = 'black') +
    geom_point(data = df_clusters, aes(x = x, y = y, color = cluster)) +
    geom_point(data = df_centers, aes(x = x, y = y), color = 'black', size = 3) +
    geom_ellipse(aes(x0 = 0, y0 = 0, a = 2.5, b = 2.5, angle = 0), color = 'red', linetype = 2, fill = 'red', alpha = 0.2) +
    geom_ellipse(aes(x0 = mu2[1], y0 = mu2[2], a = 2.5, b = 2.5, angle = 0), color = 'green', linetype = 2, fill = 'green', alpha = 0.2) +
    geom_polygon(data = df_rect, aes(x = x, y = y, group = id), fill = "blue", alpha = 0.2) +
    xlim(xlim) + ylim(ylim) + coord_fixed() +
    theme_minimal() +
    theme(panel.grid = element_blank(), legend.position = "none",
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank())
  p
}
