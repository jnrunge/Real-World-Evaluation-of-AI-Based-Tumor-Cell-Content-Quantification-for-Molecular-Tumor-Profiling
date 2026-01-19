################################################################################
# Script: 05-clustering.R
# Purpose: Perform UMAP-DBSCAN clustering to identify distinct discrepancy patterns
# Author: Jan-Niklas Runge
# 
# Description:
#   - Reduces model variable space via UMAP (Uniform Manifold Approximation)
#   - Applies DBSCAN clustering to identify sample groups with similar characteristics
#   - Performs grid search over UMAP/DBSCAN parameters to find optimal clustering
#   - Generates comprehensive visualizations: cluster embeddings, boxplots, heatmaps
#   - Analyzes term contributions and variable distributions within clusters
# 
# Outputs:
#   - output/clustering/[outcome]_*/: Parameter-specific clustering results
#   - output/clustering/best_*/: Best clustering result (highest F-ratio)
#   - *_umap_embedding.pdf: UMAP scatter plots
#   - *_dbscan_umap_clusters.pdf: Colored cluster assignments
#   - *_boxplot_discrepancy_by_cluster.pdf: Discrepancy distributions per cluster
#   - *_cluster_variable_heatmap.pdf: Variable summaries per cluster
#   - *_per_term_contributions_*.pdf: Model term contribution analyses
#   - *_parameter_grid_summary.csv: Performance metrics for all parameter combinations
#   - *_clustering_ranking.csv: Ranked clustering results by separation metric
# 
################################################################################

# Configuration ----------------------------------------------------------------
# Plot font sizes
base_size <- 14
base_size_poster <- 22

# UMAP parameters (default values; grid search explores ranges)
DEFAULT_N_NEIGHBORS <- 15  # Number of nearest neighbors for manifold approximation
DEFAULT_MIN_DIST <- 0.1    # Minimum distance between embedded points

# DBSCAN parameters (default values; grid search explores ranges)
DEFAULT_EPS <- 0.6         # Maximum distance between points in same cluster
DEFAULT_MIN_PTS_FACTOR <- 0.05  # MinPts as fraction of dataset size

# Grid search parameter ranges
EPS_VALUES <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
N_NEIGHBORS_VALUES <- c(5, 10, 15, 20)
MIN_DIST_VALUES <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

# Clustering evaluation
# F-ratio metric: between-cluster variance / within-cluster variance of tcc discrepancy
# Higher F-ratio indicates better-separated clusters

# Load models ------------------------------------------------------------------
classic_no_forced_interactions <- readRDS(file.path(project_dir, "output/processed_data/best_models_classic_no_forced_interactions.rds"))

# Helper functions -------------------------------------------------------------

#' Prepare model variables matrix
#' 
#' @param response_var Response variable name
#' @param data_df_renamed Data frame with renamed variables
#' @param model_list List of models
#' @return List with X (model matrix) and dummified_vars (factor variable names)
#' @details
#'   - Extracts variables used in model
#'   - Converts ordered factors to numeric
#'   - Creates full dummy encoding for unordered factors
prepare_model_variables <- function(response_var, data_df_renamed, model_list) {
    model_vars <- model_list[[response_var]]$model %>%
        formula() %>%
        as.character() %>%
        .[3] %>%
        strsplit(split = " + ", fixed = TRUE) %>%
        unlist() %>%
        trimws()
    
    X <- data_df_renamed %>%
        dplyr::select(all_of(model_vars[!grepl(":", model_vars, fixed = TRUE)])) %>%
        mutate(
            # characters -> factors
            across(where(is.character), as.factor),
            # ordered factors -> numeric codes
            across(where(is.ordered), ~ as.numeric(.x)),
            # logicals -> numeric (0/1) to ensure numeric matrix downstream
            across(where(is.logical), ~ as.numeric(.x))
        )

    # Identify unordered factors to be dummified
    fac_cols <- names(X)[sapply(X, is.factor)]
    dummified_vars <- if (length(fac_cols) > 0) fac_cols else NULL

    # Create full dummy variables for all remaining (unordered) factors: one column per level
    if (!is.null(dummified_vars)) {
        contr_list <- lapply(X[dummified_vars], function(x) contrasts(x, contrasts = FALSE))
        X <- as.data.frame(
            model.matrix(~ . - 1, data = X, contrasts.arg = contr_list),
            stringsAsFactors = FALSE
        )
    } else {
        # Ensure pure numeric frame if no factors present
        X <- X %>% mutate(across(where(is.integer), as.numeric))
    }

    return(list(X = X, dummified_vars = dummified_vars))
}

#' Perform UMAP transformation
#' 
#' @param X Model variables matrix
#' @param n_neighbors Number of neighbors for UMAP
#' @param min_dist Minimum distance for UMAP
#' @param seed Random seed
#' @return UMAP result as data frame with UMAP1 and UMAP2 columns
#' @details Applies scaling before UMAP to ensure equal variable influence
perform_umap <- function(X, n_neighbors = 15, min_dist = 0.1, seed = 123) {
    set.seed(seed)
    umap_result <- umap(scale(X), n_neighbors = n_neighbors, min_dist = min_dist, metric = "euclidean")
    umap_df <- as.data.frame(umap_result)
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    return(umap_df)
}

#' Plot UMAP embedding
#' 
#' @param umap_df UMAP coordinates
#' @param output_dir Output directory
#' @param response_var Response variable name
#' @param base_size Base font size for plots
plot_umap_embedding <- function(umap_df, output_dir, response_var, base_size = 14) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
        geom_point(alpha = 0.7, size = 2) +
        labs(title = "UMAP of Model Variables (scaled)") +
        theme_minimal(base_size = base_size)
    
    ggsave(file.path(output_dir, paste0(response_var, "_umap_embedding.pdf")), 
           plot = p, width = 8, height = 6, dpi = 300)
    print(p)
}

#' Determine DBSCAN parameters and perform clustering
#' 
#' @param umap_df UMAP coordinates
#' @param eps Epsilon parameter (neighborhood radius)
#' @param minPts Minimum points parameter (NULL for auto: 2 * n_dimensions)
#' @param output_dir Output directory
#' @return DBSCAN result object with cluster assignments
#' @details 
#'   - Plots k-NN distance graph to help select eps
#'   - Cluster 0 represents noise points (not assigned to any cluster)
perform_dbscan_clustering <- function(umap_df, eps = 0.6, minPts = NULL, output_dir) {
    # Plot k-NN distance for parameter selection
    kNNdist <- kNNdistplot(umap_df, k = 1:10)
    abline(h = eps, lty = 2)
    
    if (is.null(minPts)) {
        minPts <- 2 * ncol(umap_df)
    }
    db <- dbscan(umap_df, eps = eps, minPts = minPts)
    print(db)
    
    return(db)
}

#' Visualize DBSCAN clusters on UMAP
#' 
#' @param umap_df UMAP coordinates
#' @param db_clusters DBSCAN cluster assignments
#' @param output_dir Output directory
#' @param response_var Response variable name
#' @param base_size Base font size for plots
#' @param base_size_poster Base font size for poster plots
#' @return umap_df with cluster column added
#' @details Noise cluster (0) is always black; others use viridis palette
plot_dbscan_clusters <- function(umap_df, db_clusters, output_dir, response_var, 
                                 base_size = 14, base_size_poster = 22) {
    umap_df$cluster <- factor(db_clusters)
    
    # Ensure noise cluster (0) is black
    cluster_levels <- levels(umap_df$cluster)
    if ("0" %in% cluster_levels) {
        non_noise <- setdiff(cluster_levels, "0")
        vir_cols <- if (length(non_noise) > 0) viridis::viridis(length(non_noise), option = "C", begin = 0.1, end = 0.9, direction = 1) else character(0)
        color_values <- setNames(c("black", vir_cols), c("0", non_noise))
        legend_breaks <- c("0", non_noise)
    } else {
        color_values <- setNames(viridis::viridis(length(cluster_levels), option = "C", begin = 0.1, end = 0.9, direction = 1), cluster_levels)
        legend_breaks <- cluster_levels
    }

    p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_manual(values = color_values, breaks = legend_breaks, drop = FALSE) +
        labs(title = NULL, color = "Cluster") +
        theme_bw(base_size = base_size) +
        theme(legend.position = "none")
    
    print(p)
    ggsave(file.path(output_dir, paste0(response_var, "_dbscan_umap_clusters.pdf")), 
           width = 8, height = 6, dpi = 300)
    
    # Create poster version
    poster_dir <- file.path(output_dir, "poster")
    dir.create(poster_dir, showWarnings = FALSE, recursive = TRUE)
    
    p_poster <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_manual(values = color_values, breaks = legend_breaks, drop = FALSE) +
        labs(title = "UMAP of Model Variables (scaled)",
             color = "Cluster") +
        theme_bw(base_size = base_size_poster)
    
    ggsave(file.path(poster_dir, paste0(response_var, "_dbscan_umap_clusters.pdf")), 
           plot = p_poster, width = 8, height = 6, dpi = 300)
    
    return(umap_df)
}

#' Compute cluster statistics and predictions
#' 
#' @param db_clusters DBSCAN cluster assignments
#' @param response_var Response variable name
#' @param data_df Original data frame
#' @param data_df_renamed Data frame with renamed variables
#' @param model_list List of models
#' @return List with df_all (data + predictions + clusters) and cluster_counts
compute_cluster_statistics <- function(db_clusters, response_var, data_df, 
                                      data_df_renamed, model_list) {
    cluster_counts <- tibble(cluster = db_clusters) %>%
        count(cluster) %>%
        mutate(pct = n / sum(n) * 100)
    
    preds_all <- predict(
        model_list[[response_var]]$model,
        newdata = data_df_renamed
    )
    
    df_all <- data_df %>%
        mutate(
            pred = preds_all,
            res = abs(.data[[response_var]] - pred),
            cluster = db_clusters
        )
    
    return(list(df_all = df_all, cluster_counts = cluster_counts))
}

#' Create boxplot of discrepancy by cluster
#' 
#' @param df_all Data frame with clusters and predictions
#' @param cluster_counts Cluster size counts
#' @param response_var Response variable name
#' @param output_dir Output directory
#' @param base_size Base font size for plots
#' @param base_size_poster Base font size for poster plots
#' @return df_all with cluster_label column added
#' @details
#'   - Excludes noise cluster (0) from visualization
#'   - White diamonds show mean predicted value per cluster
#'   - Cluster labels include size as percentage of total data
plot_cluster_boxplots <- function(df_all, cluster_counts, response_var, output_dir, 
                                 base_size = 14, base_size_poster = 22) {
    samples_cluster_labels <- cluster_counts %>%
        mutate(cluster_label = paste0("C", cluster, "\n", round(pct, 1), "%"))
    
    df_all <- df_all %>%
        left_join(samples_cluster_labels %>% 
                     dplyr::select(cluster, cluster_label), by = "cluster")
    
    # Build color mapping: noise cluster (0) -> black, others -> viridis
    # Ensure ordering and presence of 0 label if present
    legend_breaks <- samples_cluster_labels$cluster_label
    # Identify label for cluster 0
    zero_label <- samples_cluster_labels %>%
        filter(cluster == 0) %>%
        pull(cluster_label)
    # Non-noise labels
    nn_labels <- setdiff(legend_breaks, zero_label)
    vir_cols <- if (length(nn_labels) > 0) {
        viridis::viridis(length(nn_labels), option = "C", begin = 0.1, end = 0.9, direction = 1)
    } else character(0)
    color_values <- if (length(zero_label) > 0) {
        setNames(c("black", vir_cols), c(zero_label, nn_labels))
    } else {
        setNames(vir_cols, nn_labels)
    }
    
    mean_pred_by_cluster <- df_all %>%
        dplyr::group_by(cluster_label) %>%
        summarise(mean_pred = mean(pred, na.rm = TRUE))
    
    p_summary <- ggplot(df_all %>% filter(cluster != 0), aes(x = factor(cluster_label, 
                                               levels = samples_cluster_labels$cluster_label), 
                                   y = .data[[response_var]], 
                                   fill = cluster_label)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
        geom_point(
            data = mean_pred_by_cluster %>% filter(!startsWith(cluster_label,"C0")),
            aes(x = cluster_label, y = mean_pred),
            color = "white", size = 4, shape = 18, inherit.aes = FALSE
        ) +
        labs(x = "Cluster (size % of data)",
             y = paste(response_var %>% get_pretty_name())) +
        theme_bw(base_size = base_size) +
        theme(legend.position = "none") +
        scale_fill_manual(values = color_values, breaks = legend_breaks, drop = FALSE)
    
    ggsave(file.path(output_dir, paste0(response_var, "_boxplot_discrepancy_by_cluster.pdf")), 
           plot = p_summary, width = 8, height = 6, dpi = 300)
    
    # Create poster version
    poster_dir <- file.path(output_dir, "poster")
    dir.create(poster_dir, showWarnings = FALSE, recursive = TRUE)
    
    p_summary_poster <- ggplot(df_all %>% filter(cluster != 0), aes(x = factor(cluster_label, 
                                               levels = samples_cluster_labels$cluster_label), 
                                   y = .data[[response_var]], 
                                   fill = cluster_label)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
        geom_point(
            data = mean_pred_by_cluster %>% filter(!startsWith(cluster_label,"C0")),
            aes(x = cluster_label, y = mean_pred),
            color = "white", size = 4, shape = 18, inherit.aes = FALSE
        ) +
        labs(x = "Cluster (size % of data)",
             y = paste(response_var %>% get_pretty_name())) +
        theme_bw(base_size = base_size_poster) +
        theme(legend.position = "none") +
        scale_fill_manual(values = color_values, breaks = legend_breaks, drop = FALSE)
    
    ggsave(file.path(poster_dir, paste0(response_var, "_boxplot_discrepancy_by_cluster.pdf")), 
           plot = p_summary_poster, width = 8, height = 6, dpi = 300)
    
    return(df_all)
}

#' Export extreme discrepancy counts by cluster
#' 
#' @param df_all Data frame with clusters
#' @param response_var Response variable name
#' @param output_dir Output directory
#' @details Counts samples with |discrepancy| > 20 per cluster
export_extreme_discrepancies <- function(df_all, response_var, output_dir) {
    df_all %>%
        group_by(cluster) %>%
        summarise(
            n_high = sum(.data[[response_var]] > 20, na.rm = TRUE),
            n_low = sum(.data[[response_var]] < -20, na.rm = TRUE),
            total = n(),
            pct_high = n_high / total * 100,
            pct_low = n_low / total * 100
        ) %>%
        write_excel_csv2(file.path(output_dir, paste0(response_var, "_cluster_extreme_discrepancy_counts.csv")))
}

#' Analyze term contributions per cluster
#' 
#' @param df_all Data frame with clusters
#' @param response_var Response variable name
#' @param model_list List of models
#' @param output_dir Output directory
#' @param base_size Base font size for plots
#' @details
#'   - Plots contribution of each model term to predicted value
#'   - Shows representative sample per cluster and all cluster members
analyze_term_contributions <- function(df_all, response_var, model_list, output_dir, 
                                      base_size = 14) {
    # Exclude noise cluster (0) from contribution plots
    df_all_nn <- df_all %>% filter(cluster != 0)
    if (nrow(df_all_nn) == 0 || length(unique(df_all_nn$cluster)) == 0) {
        message("No non-noise clusters to analyze for term contributions. Skipping plots.")
        return(invisible(NULL))
    }
    
    # Select representative samples
    cluster_counts <- df_all_nn %>%
        count(cluster) %>%
        mutate(pct = n / sum(n) * 100)
    
    samples <- df_all_nn %>%
        group_by(cluster) %>%
        slice_min(order_by = res, n = 1) %>%
        ungroup() %>%
        left_join(cluster_counts, by = "cluster") %>%
        mutate(
            sample_id = cluster,
            sample_label = paste0("C", cluster, "\n", round(pct, 1), "%")
        )
    
    names(samples) <- make.names(names(samples), unique = TRUE)
    
    # Get per-term contributions for samples
    contr_mat <- predict(
        model_list[[response_var]]$model,
        newdata = samples,
        type = "terms"
    )
    
    contr_list <- as_tibble(contr_mat) %>%
        mutate(sample_id = samples$sample_id) %>%
        pivot_longer(cols = -sample_id, names_to = "term", values_to = "contribution")
    
    p_contrib <- ggplot(contr_list, aes(x = term, y = contribution, 
                                        fill = contribution > 0)) +
        geom_col(show.legend = FALSE) +
        facet_wrap(~sample_id, scales = "free_x") +
        scale_x_discrete(labels = function(x) unlist(lapply(x, get_pretty_name))) +
        coord_flip() +
        labs(title = "Per-Term Contributions to the Linear Predictor") +
        theme_bw(base_size = base_size)
    
    ggsave(file.path(output_dir, paste0(response_var, "_per_term_contributions_per_sample.pdf")), 
           plot = p_contrib, width = 14, height = 6, dpi = 300)
    
    # Contributions for all cluster members
    df_all_contrib <- df_all_nn
    names(df_all_contrib) <- make.names(names(df_all_contrib), unique = TRUE)
     
    contr_mat_all <- predict(
        model_list[[response_var]]$model,
        newdata = df_all_contrib,
        type = "terms"
    )
    
    contr_list_all <- as_tibble(contr_mat_all) %>%
        mutate(cluster = df_all_contrib$cluster) %>%
        pivot_longer(cols = -cluster, names_to = "term", values_to = "contribution")
    
    contr_summary <- contr_list_all %>%
        group_by(cluster, term) %>%
        summarise(
            mean_contribution = mean(contribution, na.rm = TRUE),
            median_contribution = median(contribution, na.rm = TRUE),
            min_contribution = min(contribution, na.rm = TRUE),
            max_contribution = max(contribution, na.rm = TRUE),
            .groups = "drop"
        )
    
    p_contrib_all <- ggplot(contr_summary, aes(x = term, y = mean_contribution, 
                                               fill = cluster)) +
        geom_col(position = "dodge", show.legend = TRUE) +
        geom_errorbar(
            aes(ymin = min_contribution, ymax = max_contribution),
            position = position_dodge(width = 0.9), width = 0.2
        ) +
        facet_wrap(~cluster, scales = "free_x") +
        scale_x_discrete(labels = function(x) unlist(lapply(x, get_pretty_name))) +
        coord_flip() +
        labs(title = "Per-Term Contributions to the Linear Predictor (All Cluster Members)",
             y = "Contribution", x = "Term") +
        theme_bw(base_size = base_size)
    
    ggsave(file.path(output_dir, paste0(response_var, "_per_term_contributions_all_cluster_members.pdf")), 
           plot = p_contrib_all, width = 14, height = 8, dpi = 300)
}

#' Plot variable distributions by cluster
#' 
#' @param data_df_renamed Data frame with renamed variables
#' @param db_clusters DBSCAN cluster assignments
#' @param response_var Response variable name
#' @param model_list List of models
#' @param X Model variables matrix
#' @param output_dir Output directory
#' @param base_size Base font size for plots
#' @details
#'   - Density plots for numeric variables
#'   - Bar plots for categorical variables
#'   - Excludes noise cluster (0)
plot_variable_distributions <- function(data_df_renamed, db_clusters, response_var, 
                                       model_list, X, output_dir, base_size = 14) {
    data_df_renamed_only_model_vars <- model_list[[response_var]]$model %>%
        formula() %>%
        as.character()
    
    data_df_renamed_only_model_vars <- data_df_renamed_only_model_vars[3] %>%
        strsplit(split = " + ", fixed = TRUE) %>%
        unlist() %>%
        str_replace_all("\n", "") %>%
        str_replace_all(" ", "")
    
    clustered_data <- dplyr::select(data_df_renamed, 
                                   all_of(data_df_renamed_only_model_vars[
                                       !grepl(":", data_df_renamed_only_model_vars, fixed = TRUE)])) %>%
        mutate(cluster = factor(db_clusters))
    

    cluster_0 <- which(clustered_data$cluster == "0")

    # Exclude noise cluster (0)
    clustered_data <- clustered_data %>%
        filter(cluster != "0") %>%
        droplevels()
    
    missing_cols <- setdiff(names(X), names(clustered_data))
    if (length(missing_cols) > 0) {
        clustered_data <- bind_cols(clustered_data, X[-cluster_0, missing_cols])
    }
    
    # Numeric variables
    if (response_var %in% names(clustered_data)) {
        df_num <- clustered_data %>%
            dplyr::select(cluster, where(is.numeric), -!!response_var) %>%
            pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
    } else {
        df_num <- clustered_data %>%
            dplyr::select(cluster, where(is.numeric)) %>%
            pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
    }
    
    p_num <- ggplot(df_num, aes(x = value, fill = cluster)) +
        geom_density(alpha = 0.5) +
        facet_wrap(~variable, scales = "free", ncol = 3,
                  labeller = labeller(variable = as_labeller(
                      function(x) unlist(lapply(x, get_pretty_name))))) +
        theme_bw(base_size = base_size) +
        labs(title = "Distributions of Numeric Model Variables by Cluster",
             x = NULL, y = "Density", fill = "Cluster")
    
    ggsave(file.path(output_dir, paste0(response_var, "_numeric_variable_distributions_by_cluster.pdf")), 
           plot = p_num, width = 13, height = 8)
    
    # Categorical variables
    fac_cols <- names(clustered_data)[
        sapply(clustered_data, is.factor) & names(clustered_data) != "cluster"
    ]
    
    if (length(fac_cols) > 0) {
        df_cat <- clustered_data %>%
            dplyr::select(cluster, all_of(fac_cols)) %>%
            mutate(across(-cluster, as.character)) %>%
            pivot_longer(-cluster, names_to = "variable", values_to = "value")
        
        if (nrow(df_cat) > 0) {
            p_cat <- ggplot(df_cat, aes(x = value, fill = cluster)) +
                geom_bar(position = "dodge") +
                facet_wrap(~variable, scales = "free", ncol = 3,
                          labeller = labeller(variable = as_labeller(
                              function(x) unlist(lapply(x, get_pretty_name))))) +
                theme_bw(base_size = base_size) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                labs(title = "Distributions of Categorical Model Variables by Cluster",
                     x = NULL, y = "Count", fill = "Cluster")
            
            ggsave(file.path(output_dir, paste0(response_var, "_categorical_variable_distributions_by_cluster.pdf")), 
                   plot = p_cat, width = 12, height = 8)
        }
    }
}

#' Create cluster summary heatmap
#' 
#' @param data_df_renamed Data frame with renamed variables
#' @param db_clusters DBSCAN cluster assignments
#' @param df_all Data frame with predictions
#' @param response_var Response variable name
#' @param model_list List of models
#' @param X Model variables matrix
#' @param output_dir Output directory
#' @param dummified_vars Character vector of factor variable names
#' @param reverse_colors Logical, if TRUE blue=high, red=low
#' @param base_size Base font size for plots
#' @param base_size_poster Base font size for poster plots
#' @details
#'   - Numeric variables: z-scores relative to overall mean (++, +, =, -, --)
#'   - Categorical variables: most frequent level per cluster
#'   - Response variable: mean value per cluster with gradient coloring
#'   - Alpha channel encodes within-cluster variation
create_cluster_heatmap <- function(data_df_renamed, db_clusters, df_all, 
                                  response_var, model_list, X, output_dir,
                                  dummified_vars = NULL,
                                  reverse_colors = FALSE,
                                  base_size = 14,
                                  base_size_poster = 22) {
    # ...existing code to prepare heatmap_df...
    # ...existing code to prepare clustered_data...
    data_df_renamed_only_model_vars <- model_list[[response_var]]$model %>%
        formula() %>%
        as.character()
    
    data_df_renamed_only_model_vars <- data_df_renamed_only_model_vars[3] %>%
        strsplit(split = " + ", fixed = TRUE) %>%
        unlist() %>%
        str_replace_all("\n", "") %>%
        str_replace_all(" ", "")
    
    clustered_data <- dplyr::select(data_df_renamed, 
                                   all_of(c(dummified_vars,data_df_renamed_only_model_vars[
                                       !grepl(":", data_df_renamed_only_model_vars, fixed = TRUE)]))) %>%
        mutate(cluster = factor(db_clusters))

    # which 0s
    cluster_0 <- which(clustered_data$cluster == "0")

    # Exclude noise cluster (0)
    clustered_data <- clustered_data %>%
        filter(cluster != "0") %>%
        droplevels()
    
    missing_cols <- setdiff(names(X), names(clustered_data))

    # Remove dummified variable columns from missing_cols
        if (!is.null(dummified_vars) && length(missing_cols) > 0) {
            # Get all possible dummy column names from dummified_vars
            dummy_patterns <- unlist(lapply(dummified_vars, function(var) {
                if (var %in% names(data_df_renamed)) {
                    unique_vals <- unique(as.character(data_df_renamed[[var]]))
                    paste0(var, unique_vals)
                } else {
                    character(0)
                }
            }))
            
            # Remove missing_cols that match any dummy pattern
            missing_cols <- setdiff(missing_cols, dummy_patterns)
        }

    if (length(missing_cols) > 0) {
        clustered_data <- bind_cols(clustered_data, X[-cluster_0,missing_cols])
    }
    
    # Numeric variables summary
    if (response_var %in% names(clustered_data)) {
        df_num <- clustered_data %>%
            dplyr::select(cluster, where(is.numeric), -!!response_var) %>%
            pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
    } else {
        df_num <- clustered_data %>%
            dplyr::select(cluster, where(is.numeric)) %>%
            pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
    }
    
    num_summary <- df_num %>%
        group_by(variable, cluster) %>%
        summarise(
            mean_val = mean(value, na.rm = TRUE),
            sd_val = sd(value, na.rm = TRUE),
            n = n(),
            .groups = "drop"
        ) %>%
        left_join(
            df_num %>%
                group_by(variable) %>%
                summarise(
                    overall_mean = mean(value, na.rm = TRUE),
                    overall_sd = sd(value, na.rm = TRUE),
                    .groups = "drop"
                ),
            by = "variable"
        ) %>%
        mutate(
            z = (mean_val - overall_mean) / overall_sd,
            cat = case_when(
                z > 1.5 ~ "++",
                z > 0.5 ~ "+",
                z > -0.5 ~ "=",
                z > -1.5 ~ "-",
                TRUE ~ "--"
            ),
            variation = sd_val / overall_sd
        )
    
    # Categorical variables summary
    fac_cols <- names(clustered_data)[
        sapply(clustered_data, is.factor) & names(clustered_data) != "cluster"
    ]

    if(!is.null(dummified_vars)){
    fac_cols <- c(fac_cols, dummified_vars)
    }
    
    cat_summary <- if (length(fac_cols) > 0) {
        df_cat <- clustered_data %>%
            dplyr::select(cluster, all_of(fac_cols)) %>%
            mutate(across(-cluster, as.character)) %>%
            pivot_longer(-cluster, names_to = "variable", values_to = "value")
        
        df_cat %>%
            group_by(variable, cluster, value) %>%
            summarise(n = n(), .groups = "drop") %>%
            group_by(variable, cluster) %>%
            mutate(prop = n / sum(n)) %>%
            arrange(variable, cluster, desc(prop)) %>%
            slice(1) %>%
            ungroup() %>%
            mutate(variation = 1 - prop)
    } else {
        tibble()
    }
    
    # Response variable summary (exclude noise cluster 0)
    tc_path_summary <- df_all %>%
        filter(cluster != 0) %>%
        group_by(cluster) %>%
        summarise(
            mean_val = mean(.data[[response_var]], na.rm = TRUE),
            sd_val = sd(.data[[response_var]], na.rm = TRUE),
            n = n(),
            .groups = "drop"
        ) %>%
        mutate(
            cat = case_when(
                mean_val > 20 ~ "++",
                mean_val > 10 ~ "+",
                mean_val < -10 ~ "-",
                mean_val < -20 ~ "--",
                TRUE ~ "="
            ),
            variation = sd_val / (max(abs(mean_val), 1e-6))
        )
    
    # Prepare plot data
    num_plot_df <- num_summary %>%
        mutate(
            value_type = "numeric",
            display = cat,
            fill_legend = "Mean vs\nOverall",
            alpha_val = 1 - pmin(variation, 0.7)
        ) %>%
        select(cluster, variable, display, value_type, fill_legend, alpha_val)
    
    cat_plot_df <- cat_summary %>%
        mutate(
            value_type = "categorical",
            display = value,
            fill_legend = "Most Frequent",
            alpha_val = 1 - pmin(variation, 0.7)
        ) %>%
        select(cluster, variable, display, value_type, fill_legend, alpha_val)
    
    tc_path_plot_df <- tc_path_summary %>%
        mutate(
            variable = response_var,
            value_type = response_var,
            display = cat,
            fill_legend = "Hard Cutoff",
            alpha_val = 1 - pmin(variation, 0.7)
        ) %>%
        select(cluster, display, mean_val, variable, value_type, fill_legend, alpha_val)
    
    if (!is.factor(tc_path_plot_df$cluster)) {
        tc_path_plot_df$cluster <- factor(tc_path_plot_df$cluster)
    }
    
    heatmap_df <- bind_rows(num_plot_df, cat_plot_df, tc_path_plot_df)
    
    heatmap_df$variable <- factor(
        heatmap_df$variable,
        levels = unique(c(num_summary$variable, cat_summary$variable, response_var))
    )
    
    # Create heatmap
    p <- ggplot(heatmap_df, aes(x = factor(cluster), y = variable)) +
        geom_tile(
            data = filter(heatmap_df, value_type == "numeric") %>% 
                mutate(display = factor(display, levels = c("--", "-", "=", "+", "++"))),
            aes(fill = display, alpha = alpha_val),
            color = "white", show.legend = TRUE
        ) +
        scale_fill_manual(
            name = NULL,
            values = if (reverse_colors) {
                c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b")
            } else {
                c("#b2182b", "#ef8a62", "#f7f7f7", "#67a9cf", "#2166ac")
            },
            drop = FALSE,
            guide = guide_legend(order = 1, override.aes = list(alpha = 1))
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 1), title = NULL),
               alpha = "none") +
        scale_alpha(range = c(0.8, 1), guide = "none") +
        ggnewscale::new_scale_fill() +
        geom_tile(
            data = {
                cat_data <- subset(heatmap_df, value_type == "categorical")
                if (nrow(cat_data) > 0) {
                    # Get unique levels actually present in the data
                    present_levels <- unique(cat_data$display)
                    # Define order of levels (only those present)
                    all_possible_levels <- c("No", "Yes", "Absent", "Minimal", "Moderate", "Extensive", 
                                            "Biopsy", "Cytology", "Resection")
                    ordered_present <- all_possible_levels[all_possible_levels %in% present_levels]
                    cat_data %>% mutate(display = factor(display, levels = ordered_present))
                } else {
                    cat_data
                }
            },
            aes(fill = display, alpha = alpha_val),
            color = "white"
        ) +
        scale_fill_manual(
            name = NULL,
            values = {
                # Get levels actually present in categorical data
                cat_data <- subset(heatmap_df, value_type == "categorical")
                if (nrow(cat_data) > 0) {
                    present_levels <- unique(cat_data$display)
                    
                    # Generate palette for ordinal variables
                    pal <- grDevices::colorRampPalette(
                        if (reverse_colors) c("#2166ac", "#b2182b") else c("#b2182b", "#2166ac")
                    )(4)
                    
                    # Full color mapping
                    all_colors <- c(
                        "No" = pal[1],
                        "Yes" = pal[4],
                        "Absent" = pal[1],
                        "Minimal" = pal[2],
                        "Moderate" = pal[3],
                        "Extensive" = pal[4],
                        "Biopsy" = "#1B9E77",
                        "Cytology" = "#D95F02",
                        "Resection" = "#7570B3"
                    )
                    
                    # Return only colors for present levels
                    all_colors[names(all_colors) %in% present_levels]
                } else {
                    c()
                }
            },
            guide = guide_legend(order = 2),
            na.value = "grey90",
            drop = TRUE
        ) +
        scale_alpha(range = c(0.8, 1), guide = "none") +
        ggnewscale::new_scale_fill() +
        geom_tile(
            data = subset(heatmap_df, value_type == response_var) %>% 
                mutate(variable = get_pretty_name(response_var)),
            aes(fill = mean_val, alpha = alpha_val),
            color = "white"
        ) +
        scale_fill_gradient2(
            name = NULL,
            low = if (reverse_colors) "#2166ac" else "#b2182b",
            mid = "#f7f7f7",
            high = if (reverse_colors) "#b2182b" else "#2166ac",
            midpoint = 0,
            na.value = "grey90"
        ) +
        scale_alpha(range = c(0.8, 1), guide = "none") +
        labs(title = "Cluster-wise summary of model variables",
             x = "Cluster", y = "Variable") +
        scale_y_discrete(labels = function(x) {
            sapply(x, function(lbl) {
                pretty <- get_pretty_name(lbl)
                chars <- unlist(strsplit(pretty, ""))
                out <- ""
                count <- 0
                for (i in seq_along(chars)) {
                    out <- paste0(out, chars[i])
                    count <- count + 1
                    if (count >= 6 && chars[i] == " ") {
                        out <- paste0(out, "\n")
                        count <- 0
                    }
                }
                out
            })
        }) +
        theme_minimal(base_size = base_size) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, hjust = 0),
            panel.spacing.y = unit(0.5, "lines"),
            legend.position = "bottom"
        ) +
        coord_flip() +
        ylab(NULL)
    
    ggsave(file.path(output_dir, paste0(response_var, "_cluster_variable_heatmap.pdf")), 
           width = 16, height = 7)
    
    # Create poster version
    poster_dir <- file.path(output_dir, "poster")
    dir.create(poster_dir, showWarnings = FALSE, recursive = TRUE)
    
    p_poster <- ggplot(heatmap_df, aes(x = factor(cluster), y = variable)) +
        geom_tile(
            data = filter(heatmap_df, value_type == "numeric") %>% 
                mutate(display = factor(display, levels = c("--", "-", "=", "+", "++"))),
            aes(fill = display, alpha = alpha_val),
            color = "white", show.legend = TRUE
        ) +
        scale_fill_manual(
            name = NULL,
            values = if (reverse_colors) {
                c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b")
            } else {
                c("#b2182b", "#ef8a62", "#f7f7f7", "#67a9cf", "#2166ac")
            },
            drop = FALSE,
            guide = guide_legend(order = 1, override.aes = list(alpha = 1))
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 1), title = NULL),
               alpha = "none") +
        scale_alpha(range = c(0.8, 1), guide = "none") +
        ggnewscale::new_scale_fill() +
        geom_tile(
            data = {
                cat_data <- subset(heatmap_df, value_type == "categorical")
                if (nrow(cat_data) > 0) {
                    present_levels <- unique(cat_data$display)
                    all_possible_levels <- c("No", "Yes", "Absent", "Minimal", "Moderate", "Extensive", 
                                            "Biopsy", "Cytology", "Resection")
                    ordered_present <- all_possible_levels[all_possible_levels %in% present_levels]
                    cat_data %>% mutate(display = factor(display, levels = ordered_present))
                } else {
                    cat_data
                }
            },
            aes(fill = display, alpha = alpha_val),
            color = "white"
        ) +
        scale_fill_manual(
            name = NULL,
            values = {
                cat_data <- subset(heatmap_df, value_type == "categorical")
                if (nrow(cat_data) > 0) {
                    present_levels <- unique(cat_data$display)
                    pal <- grDevices::colorRampPalette(
                        if (reverse_colors) c("#2166ac", "#b2182b") else c("#b2182b", "#2166ac")
                    )(4)
                    all_colors <- c(
                        "No" = pal[1],
                        "Yes" = pal[4],
                        "Absent" = pal[1],
                        "Minimal" = pal[2],
                        "Moderate" = pal[3],
                        "Extensive" = pal[4],
                        "Biopsy" = "#1B9E77",
                        "Cytology" = "#D95F02",
                        "Resection" = "#7570B3"
                    )
                    all_colors[names(all_colors) %in% present_levels]
                } else {
                    c()
                }
            },
            guide = guide_legend(order = 2),
            na.value = "grey90",
            drop = TRUE
        ) +
        scale_alpha(range = c(0.8, 1), guide = "none") +
        ggnewscale::new_scale_fill() +
        geom_tile(
            data = subset(heatmap_df, value_type == response_var) %>% 
                mutate(variable = get_pretty_name(response_var)),
            aes(fill = mean_val, alpha = alpha_val),
            color = "white"
        ) +
        scale_fill_gradient2(
            name = NULL,
            low = if (reverse_colors) "#2166ac" else "#b2182b",
            mid = "#f7f7f7",
            high = if (reverse_colors) "#b2182b" else "#2166ac",
            midpoint = 0,
            na.value = "grey90"
        ) +
        scale_alpha(range = c(0.8, 1), guide = "none") +
        labs(title = "Cluster-wise summary of model variables",
             x = "Cluster", y = "Variable") +
        scale_y_discrete(labels = function(x) {
            sapply(x, function(lbl) {
                pretty <- get_pretty_name(lbl)
                chars <- unlist(strsplit(pretty, ""))
                out <- ""
                count <- 0
                for (i in seq_along(chars)) {
                    out <- paste0(out, chars[i])
                    count <- count + 1
                    if (count >= 6 && chars[i] == " ") {
                        out <- paste0(out, "\n")
                        count <- 0
                    }
                }
                out
            })
        }) +
        theme_minimal(base_size = base_size_poster) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, hjust = 0),
            panel.spacing.y = unit(0.5, "lines"),
            legend.position = "bottom"
        ) +
        coord_flip() +
        ylab(NULL)
    
    ggsave(file.path(poster_dir, paste0(response_var, "_cluster_variable_heatmap.pdf")), 
           plot = p_poster, width = 20, height = 10)
}

#' Calculate cluster separation metric for TCC discrepancy
#' 
#' @param df_all Data frame with cluster assignments and TCC discrepancy
#' @param response_var Response variable name
#' @return List with separation metrics (f_ratio, between_var, within_var, etc.)
#' @details
#'   - F-ratio = between-cluster variance / within-cluster variance
#'   - Higher F-ratio indicates better separation
#'   - Used to select best clustering from grid search
calculate_cluster_separation <- function(df_all, response_var) {
    if (is.null(df_all) || !"cluster" %in% names(df_all)) {
        return(list(f_ratio = NA, between_var = NA, within_var = NA, n_clusters = NA))
    }
    
    # Calculate overall mean and variance
    overall_mean <- mean(df_all[[response_var]], na.rm = TRUE)
    total_var <- var(df_all[[response_var]], na.rm = TRUE)
    
    # Calculate between-cluster and within-cluster variance
    cluster_stats <- df_all %>%
        group_by(cluster) %>%
        summarise(
            cluster_mean = mean(.data[[response_var]], na.rm = TRUE),
            cluster_var = var(.data[[response_var]], na.rm = TRUE),
            n = n(),
            .groups = "drop"
        )
    
    # Between-cluster variance (weighted by cluster size)
    between_var <- sum(cluster_stats$n * (cluster_stats$cluster_mean - overall_mean)^2) / 
                   sum(cluster_stats$n)
    
    # Within-cluster variance (weighted average)
    within_var <- sum(cluster_stats$n * cluster_stats$cluster_var, na.rm = TRUE) / 
                  sum(cluster_stats$n)
    
    # F-ratio (higher is better separation)
    f_ratio <- between_var / (within_var + 1e-10)
    
    # Number of non-noise clusters
    n_clusters <- sum(cluster_stats$cluster != 0)
    
    return(list(
        f_ratio = f_ratio,
        between_var = between_var,
        within_var = within_var,
        total_var = total_var,
        n_clusters = n_clusters,
        cluster_means = cluster_stats$cluster_mean,
        cluster_sizes = cluster_stats$n
    ))
}

#' Find best clustering result from grid search
#' 
#' @param grid_results Results from run_clustering_grid_search
#' @param response_var Response variable name
#' @return List with ranking, best_result, best_index, best_params
#' @details
#'   - Ranks all parameter combinations by F-ratio
#'   - Renames output directory of best result to include "best_" prefix
find_best_clustering <- function(grid_results, response_var) {
    message("Evaluating ", length(grid_results$results), " clustering results...")
    
    # Calculate separation metrics for each result
    separation_metrics <- lapply(seq_along(grid_results$results), function(i) {
        result <- grid_results$results[[i]]
        
        if ("error" %in% names(result) || is.null(result$df_all)) {
            return(tibble(
                index = i,
                f_ratio = NA,
                between_var = NA,
                within_var = NA,
                n_clusters = NA,
                has_error = TRUE
            ))
        }
        
        metrics <- calculate_cluster_separation(result$df_all, response_var)
        
        tibble(
            index = i,
            f_ratio = metrics$f_ratio,
            between_var = metrics$between_var,
            within_var = metrics$within_var,
            total_var = metrics$total_var,
            n_clusters = metrics$n_clusters,
            has_error = FALSE
        )
    })
    
    # Combine with parameters
    ranking <- bind_rows(separation_metrics) %>%
        bind_cols(grid_results$parameters) %>%
        arrange(desc(f_ratio))
    
    # Get best result
    best_idx <- ranking$index[1]
    best_result <- grid_results$results[[best_idx]]
    best_params <- grid_results$parameters[best_idx, ]
    
    # Rename output directory to include "best_" prefix
    old_output_dir <- best_result$output_dir
    dir_name <- basename(old_output_dir)
    parent_dir <- dirname(old_output_dir)
    new_dir_name <- paste0("best_", dir_name)
    new_output_dir <- file.path(parent_dir, new_dir_name)
    
    # Rename the directory
    if (dir.exists(old_output_dir)) {
        # Remove new directory if it already exists
        if (dir.exists(new_output_dir)) {
            unlink(new_output_dir, recursive = TRUE)
        }
        file.rename(old_output_dir, new_output_dir)
        best_result$output_dir <- new_output_dir
        message("Renamed output directory to: ", new_output_dir)
    }
    
    message("\n========================================")
    message("BEST CLUSTERING RESULT:")
    message("Index: ", best_idx)
    message("Parameters: eps=", best_params$eps, 
            ", minPts=", best_params$minPts,
            ", n_neighbors=", best_params$n_neighbors,
            ", min_dist=", best_params$min_dist)
    message("F-ratio: ", round(ranking$f_ratio[1], 3))
    message("Number of clusters: ", ranking$n_clusters[1])
    message("Output directory: ", new_output_dir)
    message("========================================\n")
    
    # Print top 5 results
    message("Top 5 parameter combinations:")
    print(ranking %>% 
          select(index, eps, minPts, n_neighbors, min_dist, f_ratio, n_clusters) %>%
          head(5))
    
    return(list(
        ranking = ranking,
        best_result = best_result,
        best_index = best_idx,
        best_params = best_params
    ))
}

# Main workflow function -------------------------------------------------------

#' Main clustering analysis workflow
#' 
#' @param response_var Response variable name
#' @param data_df Original data frame
#' @param data_df_renamed Data frame with renamed variables
#' @param model_list List containing models
#' @param output_base_dir Base output directory
#' @param eps DBSCAN epsilon parameter
#' @param minPts DBSCAN minimum points (default: 2 * UMAP dimensions)
#' @param n_neighbors UMAP n_neighbors parameter
#' @param min_dist UMAP min_dist parameter
#' @param seed Random seed
#' @param reverse_colors Logical, reverse heatmap colors
#' @param base_size Base font size for plots
#' @param base_size_poster Base font size for poster plots
#' @return List with clustering results
run_clustering_analysis <- function(response_var, 
                                   data_df, 
                                   data_df_renamed, 
                                   model_list,
                                   output_base_dir = file.path(project_dir, "output/clustering"),
                                   eps = 0.6,
                                   minPts = NULL,
                                   n_neighbors = 15,
                                   min_dist = 0.1,
                                   seed = 123,
                                   reverse_colors = FALSE,
                                   base_size = 14,
                                   base_size_poster = 22) {
    
    # Create parameter-specific output directory
    param_dir <- sprintf("%s_eps%.2f_minPts%d_nn%d_md%.2f", 
                        response_var, eps, 
                        ifelse(is.null(minPts), 4, minPts), 
                        n_neighbors, min_dist)
    output_dir <- file.path(output_base_dir, param_dir)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    message("Starting clustering analysis for: ", response_var)
    message("Parameters: eps=", eps, ", minPts=", ifelse(is.null(minPts), "auto", minPts), 
            ", n_neighbors=", n_neighbors, ", min_dist=", min_dist)
    message("Output directory: ", output_dir)
    
    # 1. Prepare model variables
    message("Step 1: Preparing model variables...")
    vars <- prepare_model_variables(response_var, data_df_renamed, model_list)
    X <- vars$X
    dummified_vars <- vars$dummified_vars

    # 2. UMAP transformation
    message("Step 2: Performing UMAP transformation...")
    umap_df <- perform_umap(X, n_neighbors = n_neighbors, min_dist = min_dist, seed = seed)
    plot_umap_embedding(umap_df, output_dir, response_var, base_size = base_size)
    
    # 3. DBSCAN clustering
    message("Step 3: Performing DBSCAN clustering...")
    db <- perform_dbscan_clustering(umap_df, eps = eps, minPts = minPts, output_dir = output_dir)
    
    # 4. Visualize clusters
    message("Step 4: Visualizing clusters...")
    umap_df <- plot_dbscan_clusters(umap_df, db$cluster, output_dir, response_var, 
                                    base_size = base_size, base_size_poster = base_size_poster)
    
    # 5. Compute cluster statistics
    message("Step 5: Computing cluster statistics...")
    cluster_stats <- compute_cluster_statistics(db$cluster, response_var, 
                                                data_df, data_df_renamed, model_list)
    df_all <- cluster_stats$df_all
    cluster_counts <- cluster_stats$cluster_counts
    
    # 6. Plot cluster boxplots
    message("Step 6: Creating cluster boxplots...")
    df_all <- plot_cluster_boxplots(df_all, cluster_counts, response_var, output_dir,
                                    base_size = base_size, base_size_poster = base_size_poster)
    
    # 7. Export extreme discrepancies
    message("Step 7: Exporting extreme discrepancy counts...")
    export_extreme_discrepancies(df_all, response_var, output_dir)
    
    # 8. Analyze term contributions
    message("Step 8: Analyzing term contributions...")
    analyze_term_contributions(df_all, response_var, model_list, output_dir, base_size = base_size)
    
    # 9. Plot variable distributions
    message("Step 9: Plotting variable distributions...")
    plot_variable_distributions(data_df_renamed, db$cluster, response_var, 
                                model_list, X, output_dir, base_size = base_size)
    
    # 10. Create cluster heatmap
    message("Step 10: Creating cluster summary heatmap...")
    create_cluster_heatmap(data_df_renamed, db$cluster, df_all, response_var, 
                           model_list, X, output_dir, 
                           dummified_vars = dummified_vars,
                           reverse_colors = reverse_colors,
                           base_size = base_size,
                           base_size_poster = base_size_poster)

    message("Clustering analysis complete!")
    
    return(list(
        clusters = db$cluster,
        df_all = df_all,
        cluster_counts = cluster_counts,
        umap_df = umap_df,
        parameters = list(eps = eps, minPts = minPts, n_neighbors = n_neighbors, min_dist = min_dist),
        output_dir = output_dir
    ))
}

#' Run clustering analysis over parameter grid
#' 
#' @param response_var Response variable name
#' @param data_df Original data frame
#' @param data_df_renamed Data frame with renamed variables
#' @param model_list List containing models
#' @param output_base_dir Base output directory
#' @param eps_values Vector of epsilon values to try
#' @param minPts_values Vector of minPts values to try (NULL for auto)
#' @param n_neighbors_values Vector of n_neighbors values to try
#' @param min_dist_values Vector of min_dist values to try
#' @param seed Random seed
#' @param n_cores Number of cores to use (default: detectCores() - 1)
#' @param reverse_colors Logical, reverse heatmap colors
#' @param base_size Base font size for plots
#' @param base_size_poster Base font size for poster plots
#' @return List with results, parameters, and summary
#' @details
#'   - Runs clustering for all parameter combinations in parallel
#'   - Saves summary table with cluster counts and error status
run_clustering_grid_search <- function(response_var,
                                      data_df,
                                      data_df_renamed,
                                      model_list,
                                      output_base_dir = file.path(project_dir, "output/clustering"),
                                      eps_values = c(0.5, 0.6, 0.7),
                                      minPts_values = c(5, 10, 15),
                                      n_neighbors_values = c(10, 15, 20),
                                      min_dist_values = c(0.05, 0.1, 0.2),
                                      seed = 123,
                                      n_cores = NULL,
                                      reverse_colors = FALSE,
                                      base_size = 14,
                                      base_size_poster = 22) {
    
    # Create parameter grid
    param_grid <- expand.grid(
        eps = eps_values,
        minPts = minPts_values,
        n_neighbors = n_neighbors_values,
        min_dist = min_dist_values,
        stringsAsFactors = FALSE
    )
    
    # Determine number of cores
    if (is.null(n_cores)) {
        n_cores <- max(1, detectCores() - 1)
    }
    n_cores <- min(n_cores, nrow(param_grid))
    
    message("Running clustering analysis over ", nrow(param_grid), " parameter combinations...")
    message("Using ", n_cores, " cores for parallel processing")
    
    # Setup parallel backend
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Export necessary objects and functions to cluster
    clusterExport(cl, c("response_var", "data_df", "data_df_renamed", "model_list", 
                       "output_base_dir", "seed", "reverse_colors", "base_size", "base_size_poster",
                       "run_clustering_analysis", "prepare_model_variables",
                       "perform_umap", "plot_umap_embedding", "perform_dbscan_clustering",
                       "plot_dbscan_clusters", "compute_cluster_statistics",
                       "plot_cluster_boxplots", "export_extreme_discrepancies",
                       "analyze_term_contributions", "plot_variable_distributions",
                       "create_cluster_heatmap", "get_pretty_name", "pretty_names"),
                 envir = environment())
    
    # Load required packages on each worker
    clusterEvalQ(cl, {
        library(dbscan)
        library(ggplot2)
        library(viridis)
        library(factoextra)
        library(dplyr)
        library(uwot)
        library(tidyr)
        library(readr)
        library(stringr)
        library(ggnewscale)
    })
    
    # Run parallel grid search
    all_results <- tryCatch({
        foreach(i = 1:nrow(param_grid), 
                .packages = c("dplyr", "dbscan", "uwot", "ggplot2"),
                .errorhandling = "pass") %dopar% {
            
            params <- param_grid[i, ]
            
            tryCatch({
                result <- run_clustering_analysis(
                    response_var = response_var,
                    data_df = data_df,
                    data_df_renamed = data_df_renamed,
                    model_list = model_list,
                    output_base_dir = output_base_dir,
                    eps = params$eps,
                    minPts = params$minPts,
                    n_neighbors = params$n_neighbors,
                    min_dist = params$min_dist,
                    seed = seed,
                    reverse_colors = reverse_colors,
                    base_size = base_size,
                    base_size_poster = base_size_poster
                )
                result
            }, error = function(e) {
                list(error = e$message, parameters = params)
            })
        }
    }, finally = {
        stopCluster(cl)
    })
    
    # Create summary
    summary_df <- param_grid %>%
        mutate(
            n_clusters = sapply(all_results, function(r) {
                if ("clusters" %in% names(r)) length(unique(r$clusters)) else NA
            }),
            has_error = sapply(all_results, function(r) "error" %in% names(r)),
            output_dir = sapply(all_results, function(r) {
                if ("output_dir" %in% names(r)) r$output_dir else NA_character_
            })
        )
    
    write_csv(summary_df, file.path(output_base_dir, paste0(response_var, "_parameter_grid_summary.csv")))
    
    message("\n========================================")
    message("Grid search complete!")
    message("Summary saved to: ", file.path(output_base_dir, paste0(response_var, "_parameter_grid_summary.csv")))
    message("========================================")
    
    return(list(
        results = all_results,
        parameters = param_grid,
        summary = summary_df
    ))
}

# Execute analysis -------------------------------------------------------------
message("\n=== Starting clustering analysis ===")

message("Found ", length(response_vars), " response variables to analyze:")
message(paste(response_vars, collapse = ", "))

# Loop through all response variables
for (response_var in response_vars) {
    message("\n\n========================================")
    message("STARTING ANALYSIS FOR: ", response_var)
    message("========================================\n")
    
    # Run grid search with parameter ranges defined in configuration
    grid_results <- run_clustering_grid_search(
        response_var = response_var,
        data_df = data_df,
        data_df_renamed = data_df_renamed,
        model_list = classic_no_forced_interactions,
        output_base_dir = file.path(project_dir, "output/clustering"),
        eps_values = EPS_VALUES,
        minPts_values = nrow(data_df_renamed) * DEFAULT_MIN_PTS_FACTOR,
        n_neighbors_values = N_NEIGHBORS_VALUES,
        min_dist_values = MIN_DIST_VALUES,
        seed = 123,
        reverse_colors = TRUE  # Blue=high, red=low for heatmaps
    )
    
    # Save grid results for future reference
    saveRDS(grid_results, 
            file = file.path(project_dir, "output/clustering", 
                             paste0(response_var, "_clustering_grid_results.rds")))
    
    # Find and rename best clustering result
    best_clustering <- find_best_clustering(grid_results, response_var)
    
    # Save ranking of all parameter combinations
    write_csv(
        best_clustering$ranking,
        file.path(project_dir, "output/clustering",
            paste0(response_var, "_clustering_ranking.csv")
        )
    )
    
    message("\n========================================")
    message("COMPLETED ANALYSIS FOR: ", response_var)
    message("========================================\n")
}

message("\n\n========================================")
message("ALL RESPONSE VARIABLES PROCESSED")
message("Results saved to:", file.path(project_dir, "output/clustering"))
message("========================================")