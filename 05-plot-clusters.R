## cluster and plot

vars_in_model <- all_types_not_enforced_interaction$model$coefficients[-1] %>% names()

vars_in_model <- str_replace(vars_in_model, "\\.(L|C|Q)$", "") %>% unique()

set.seed(123)
# 1. Cluster on the explanatory variables
X <- all_model_variables_renamed %>%
    dplyr::select(all_of(vars_in_model)) %>%
    mutate(
        across(where(is.character), as.factor),
        across(where(is.factor), as.numeric),
        across(where(is.numeric), ~ scale(.)[, 1])
    )

# PCA for 2D reduction
pca <- prcomp(X, center = TRUE, scale. = TRUE)
X_pca <- as_tibble(pca$x[, 1:2]) %>%
    rename(PC1 = 1, PC2 = 2)

# K-means clustering
set.seed(123)
km <- kmeans(scale(X), centers = 5)

# Add cluster assignments to PCA data
X_pca <- X_pca %>%
    mutate(cluster = factor(km$cluster))

cluster_counts <- tibble(
    cluster = 1:length(km$size),
    count   = as.integer(km$size)
) %>%
    mutate(pct = count / sum(count) * 100)

# Example plot (optional)
ggplot(X_pca, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(size = 2, alpha = 0.7) +
    theme_bw() +
    labs(title = "PCA of Model Variables with K-means Clusters")

X_pca$TC_Path_minus_TC_AI <- all_model_variables_renamed$TC_Path_minus_TC_AI

ggplot(X_pca, aes(x = PC1, y = PC2, color = TC_Path_minus_TC_AI)) +
    geom_point(size = 5, alpha = 0.7, shape = 21, stroke = 0.7, fill = NA) +
    geom_point(aes(fill = TC_Path_minus_TC_AI), size = 5, alpha = 0.7, shape = 21, color = "black", stroke = 0.7) +
    theme_bw() +
    labs(
        title = "PCA of samples with TCC Discrepancy",
        color = "TC Path - TC AI"
    ) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-100, 100)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-100, 100))


# --- PCA heatmap of mean TCC discrepancy ---

# Control the grid resolution (smaller = finer, larger = coarser)
pca_grid_size <- 1 # in PCA units; increase for more merging

# Compute grid boundaries
pc1_range <- range(X_pca$PC1, na.rm = TRUE)
pc2_range <- range(X_pca$PC2, na.rm = TRUE)

# Create grid centers
pc1_seq <- seq(pc1_range[1], pc1_range[2], length.out = ceiling((pc1_range[2] - pc1_range[1]) / pca_grid_size))
pc2_seq <- seq(pc2_range[1], pc2_range[2], length.out = ceiling((pc2_range[2] - pc2_range[1]) / pca_grid_size))

# Assign each point to a grid cell
X_pca <- X_pca %>%
    mutate(
        grid_PC1 = cut(PC1, breaks = c(pc1_seq - pca_grid_size / 2, max(pc1_seq) + pca_grid_size / 2), labels = pc1_seq),
        grid_PC2 = cut(PC2, breaks = c(pc2_seq - pca_grid_size / 2, max(pc2_seq) + pca_grid_size / 2), labels = pc2_seq)
    )

# Compute mean discrepancy per grid cell
heatmap_df <- X_pca %>%
    filter(!is.na(grid_PC1) & !is.na(grid_PC2)) %>%
    group_by(grid_PC1, grid_PC2) %>%
    summarise(
        mean_discrepancy = mean(TC_Path_minus_TC_AI, na.rm = TRUE),
        n = n(),
        .groups = "drop"
    ) %>%
    mutate(
        grid_PC1 = as.numeric(as.character(grid_PC1)),
        grid_PC2 = as.numeric(as.character(grid_PC2))
    )

# Compute cluster centers and radii in PCA space
cluster_centers <- as.data.frame(km$centers %*% pca$rotation[, 1:2])
colnames(cluster_centers) <- c("PC1", "PC2")
cluster_centers$cluster <- factor(1:nrow(cluster_centers))

# For each cluster, compute the maximum distance from center to any point in that cluster (in PCA space)
# Compute cluster radii as a quantile (e.g., 80th percentile) of distances to center, for tighter circles
cluster_radii <- sapply(1:nrow(cluster_centers), function(k) {
    pts <- X_pca %>% filter(cluster == k)
    if (nrow(pts) == 0) {
        return(0)
    }
    dists <- sqrt((pts$PC1 - cluster_centers$PC1[k])^2 + (pts$PC2 - cluster_centers$PC2[k])^2)
    quantile(dists, 0.8, na.rm = TRUE) # 80th percentile distance
})
cluster_centers$radius <- cluster_radii

# Plot as a heatmap (geom_tile), overlay cluster circles
ggplot() +
    geom_tile(
        data = heatmap_df,
        aes(x = grid_PC1, y = grid_PC2, fill = mean_discrepancy),
        color = NA
    ) +
    geom_circle(
        data = cluster_centers,
        aes(x0 = PC1, y0 = PC2, r = radius, color = cluster),
        fill = NA,
        size = 0.5,
        inherit.aes = FALSE
    ) +
    geom_point(
        data = X_pca,
        aes(x = PC1, y = PC2),
        size = 2, shape = 1, alpha = 1
    ) +
    scale_fill_gradient2(
        low = "red", mid = "white", high = "blue",
        midpoint = 0, limits = c(-100, 100), na.value = "grey80"
    ) +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.background = element_blank()
    ) +
    labs(
        title = "PCA Heatmap of Mean TCC Discrepancy",
        fill = "Mean TC Path - TC AI"
    )



ggplot(X_pca, aes(cluster, TC_Path_minus_TC_AI)) +
    geom_boxplot() +
    theme_bw() +
    labs(
        title = "TCC Discrepancy by Cluster",
        x = "Cluster",
        y = "TC Path - TC AI"
    ) +
    scale_x_discrete(labels = paste0("Cluster ", 1:5))

# Show distribution of each variable in X by cluster
X_long <- X %>%
    mutate(cluster = X_pca$cluster) %>%
    pivot_longer(
        cols = -cluster,
        names_to = "variable",
        values_to = "value"
    )
# Apply get_pretty_name to all variable names for labelling
pretty_names_vec <- setNames(
    vapply(unique(X_long$variable), get_pretty_name, character(1)),
    unique(X_long$variable)
)

ggplot(X_long, aes(x = cluster, y = value, fill = cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    facet_wrap(
        ~variable,
        scales = "free_y",
        ncol = 4,
        labeller = as_labeller(pretty_names_vec)
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
        title = "Distribution of Model Variables by Cluster",
        x = "Cluster",
        y = "Scaled Value"
    )
