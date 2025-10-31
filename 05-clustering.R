library(dbscan)
library(ggplot2)
library(viridis) # for colorblind-friendly palettes
library(factoextra)
library(dplyr)
library(uwot)


## 1. pick 5 “far-apart” samples by k-means on the predictors
set.seed(123)
# 1. Cluster on the explanatory variables
# Only include variables present in the model formula
model_vars <- classic_no_forced_interactions$model %>%
    formula() %>%
    as.character() %>%
    .[3] %>%
    strsplit(split = " + ", fixed = TRUE) %>%
    unlist() %>%
    trimws()

# 2. compute predictions and residuals



X <- data_df_renamed %>%
    dplyr::select(all_of(model_vars[!grepl(":",model_vars,fixed=TRUE)])) %>% 
    mutate(
        across(where(is.character), as.factor),
        across(where(is.factor), as.numeric)
    )


## dbscan

# decide on epsiolon an minpts

X2 <- bind_cols(X, dplyr::select(data_df_renamed, TCC_Patho_minus_TCC_AI))

# UMAP transform first

set.seed(123)
umap_result <- umap(scale(X), n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Plot UMAP embedding colored by cluster
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(alpha = 0.7, size = 2) +
    labs(title = "UMAP of Model Variables (scaled)") +
    theme_minimal()

kNNdist <- kNNdistplot(umap_df, k = 1:10)
abline(h = 0.5, lty = 2)

minpts <- 2 * ncol(umap_df)

db <- dbscan(umap_df, eps = 0.5, minPts = minpts)
print(db)

umap_df$cluster <- factor(db$cluster)



ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_viridis_d(option = "C", begin = 0.1, end = 0.9, direction = 1) +
    labs(title = "UMAP of Model Variables (scaled) colored by DBSCAN cluster") +
    theme_minimal(14) +
    labs(color = "Cluster")


dir.create("discrepancies/2025-10-Data-Version/clustering", showWarnings = FALSE)
ggsave("discrepancies/2025-10-Data-Version/clustering/dbscan_umap_clusters.png", width = 8, height = 6, dpi = 300)















# 1a. compute cluster sizes (as % of all rows)
cluster_counts <- tibble(cluster = db$cluster) %>%
    count(cluster) %>%
    mutate(pct = n / sum(n) * 100)

preds_all <- predict(
    classic_no_forced_interactions$model,
    newdata = data_df_renamed
)

df_all <- data_df %>%
    mutate(
        pred    = preds_all,
        res     = abs(TCC_Patho_minus_TCC_AI - pred),
        cluster = db$cluster
    )

# 3. pick one “well‐explained” sample per cluster
samples <- df_all %>%
    group_by(cluster) %>%
    slice_min(order_by = res, n = 1) %>%
    ungroup() %>%
    # join cluster size pct and assign sample_id
    left_join(cluster_counts, by = "cluster") %>%
    mutate(
        sample_id    = cluster,
        sample_label = paste0("C", cluster, "\n", round(pct, 1), "%")
    )

# 4. ensure syntactic names
names(samples) <- make.names(names(samples), unique = TRUE)

# 5. get predictions + 95% CI for these samples
# Boxplots of observed TC_Path_minus_TC_AI values for the four clusters

# Add cluster labels with size percentage
samples_cluster_labels <- cluster_counts %>%
    mutate(cluster_label = paste0("C", cluster, "\n", round(pct, 1), "%"))

df_all <- df_all %>%
    left_join(samples_cluster_labels %>% dplyr::select(cluster, cluster_label), by = "cluster")

# Compute mean predicted value per cluster
mean_pred_by_cluster <- df_all %>%
    dplyr::group_by(cluster_label) %>%
    summarise(mean_pred = mean(pred, na.rm = TRUE))

# Plot boxplots of observed values by cluster, add mean predicted value as a point
p_summary <- ggplot(df_all, aes(x = factor(cluster_label, levels = samples_cluster_labels$cluster_label), y = TCC_Patho_minus_TCC_AI, fill = cluster_label)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    geom_point(
        data = mean_pred_by_cluster,
        aes(x = cluster_label, y = mean_pred),
        color = "white", size = 4, shape = 18, inherit.aes = FALSE
    ) +
    labs(
        x = "Cluster (size % of data)",
        y = "TCC Path-AI discrepancy",
        title = NULL,
        subtitle = NULL
    ) +
    theme_bw(14) +
    theme(legend.position = "none") +
    scale_fill_viridis_d(option = "C", begin = 0.1, end = 0.9, direction = 1)

ggsave("discrepancies/2025-10-Data-Version/clustering/boxplot_TCC_discrepancy_by_cluster.png", plot = p_summary, width = 8, height = 6, dpi = 300)











df_all %>%
    group_by(cluster) %>%
    summarise(
        n_high = sum(TCC_Patho_minus_TCC_AI > 20, na.rm = TRUE),
        n_low = sum(TCC_Patho_minus_TCC_AI < -20, na.rm = TRUE),
        total = n(),
        pct_high = n_high / total * 100,
        pct_low = n_low / total * 100
    ) %>%
    print()

## 4. decompose each prediction into covariate contributions
# get each term’s contribution via predict(type="terms")
contr_mat <- predict(
    classic_no_forced_interactions$model,
    newdata = samples,
    type = "terms"
)
contr_list <- as_tibble(contr_mat) %>%
    mutate(sample_id = samples$sample_id) %>%
    pivot_longer(
        cols = -sample_id,
        names_to = "term",
        values_to = "contribution"
    )

# now contr_list has one row per sample_id×term
head(contr_list)

p_contrib <- ggplot(contr_list, aes(x = term, y = contribution, fill = contribution > 0)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~sample_id, scales = "free_x") +
    scale_x_discrete(
        labels = function(x) {
            unlist(lapply(x, get_pretty_name))
        }
    ) +
    coord_flip() +
    labs(title = "Per‐Term Contributions to the Linear Predictor") +
    theme_bw(base_size = 11)

## 5. print plots and a one‐liner for each

ggsave("discrepancies/2025-10-Data-Version/clustering/per_term_contributions_per_sample.png", plot = p_contrib, width = 8, height = 6, dpi = 300)

df_all_contrib <- df_all
names(df_all_contrib) <- make.names(names(df_all_contrib), unique = TRUE)


# Calculate contributions for all members of each cluster
contr_mat_all <- predict(
    classic_no_forced_interactions$model,
    newdata = df_all_contrib, # Use all data instead of just samples
    type = "terms"
)

# Convert contributions to a tidy format
contr_list_all <- as_tibble(contr_mat_all) %>%
    mutate(cluster = df_all_contrib$cluster) %>%
    pivot_longer(
        cols = -cluster,
        names_to = "term",
        values_to = "contribution"
    )

# Summarize contributions for each term within each cluster
contr_summary <- contr_list_all %>%
    group_by(cluster, term) %>%
    summarise(
        mean_contribution = mean(contribution, na.rm = TRUE),
        median_contribution = median(contribution, na.rm = TRUE),
        min_contribution = min(contribution, na.rm = TRUE),
        max_contribution = max(contribution, na.rm = TRUE),
        .groups = "drop"
    )

# Visualize the range of contributions for each term within each cluster
p_contrib_all <- ggplot(contr_summary, aes(x = term, y = mean_contribution, fill = cluster)) +
    geom_col(position = "dodge", show.legend = TRUE) +
    geom_errorbar(
        aes(
            ymin = min_contribution,
            ymax = max_contribution
        ),
        position = position_dodge(width = 0.9),
        width = 0.2
    ) +
    facet_wrap(~cluster, scales = "free_x") +
    scale_x_discrete(
        labels = function(x) {
            unlist(lapply(x, get_pretty_name))
        }
    ) +
    coord_flip() +
    labs(
        title = "Per-Term Contributions to the Linear Predictor (All Cluster Members)",
        y = "Contribution",
        x = "Term"
    ) +
    theme_minimal(base_size = 11)

ggsave("discrepancies/2025-10-Data-Version/clustering/per_term_contributions_all_cluster_members.png", plot = p_contrib_all, width = 10, height = 8, dpi = 300)


















data_df_renamed_only_model_vars <- classic_no_forced_interactions$model %>%
    formula() %>%
    as.character()

data_df_renamed_only_model_vars <- data_df_renamed_only_model_vars[3] %>%
    strsplit(split = " + ", fixed = TRUE) %>%
    unlist()

data_df_renamed_only_model_vars <- data_df_renamed_only_model_vars %>%
    str_replace_all("\n", "") %>%
    str_replace_all(" ", "")

clustered_data <- dplyr::select(data_df_renamed, all_of(data_df_renamed_only_model_vars[!grepl(":",data_df_renamed_only_model_vars,fixed=TRUE)])) %>%
    mutate(cluster = factor(db$cluster))

# Ensure all columns in X are present in clustered_data
missing_cols <- setdiff(names(X), names(clustered_data))
if (length(missing_cols) > 0) {
    clustered_data <- bind_cols(clustered_data, X[missing_cols])
}

# 1) Numeric variables: density by cluster
if (response_var %in% names(clustered_data)) {
    df_num <- clustered_data %>%
        dplyr::select(cluster, where(is.numeric), -!!response_var) %>%
        pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
} else {
    df_num <- clustered_data %>%
        dplyr::select(cluster, where(is.numeric)) %>%
        pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
}

# Special handling for TC_Path_minus_TC_AI if present
has_tc_path <- response_var %in% names(clustered_data)
if (has_tc_path) {
    df_tc_path <- clustered_data %>%
        dplyr::select(cluster, !!response_var) %>%
        pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
}

p_num <- ggplot(df_num, aes(x = value, fill = cluster)) +
    geom_density(alpha = 0.5) +
    facet_wrap(
        ~variable,
        scales = "free",
        ncol = 3,
        labeller = labeller(variable = as_labeller(function(x) {
            unlist(lapply(x, get_pretty_name))
        }))
    ) +
    theme_bw() +
    labs(
        title = "Distributions of Numeric Model Variables by Cluster",
        x = NULL, y = "Density", fill = "Cluster"
    )


ggsave("discrepancies/2025-10-Data-Version/clustering/numeric_variable_distributions_by_cluster.pdf", plot = p_num, width = 12, height = 8)



# 2) Categorical variables: bar chart by cluster
fac_cols <- names(clustered_data)[
    sapply(clustered_data, is.factor) & names(clustered_data) != "cluster"
]
if (length(fac_cols) > 0) {
    df_cat <- clustered_data %>%
        dplyr::select(cluster, all_of(fac_cols)) %>%
        mutate(across(-cluster, as.character)) %>% # convert all factor columns to character
        pivot_longer(-cluster, names_to = "variable", values_to = "value")

    if (nrow(df_cat) > 0) {
        p_cat <- ggplot(df_cat, aes(x = value, fill = cluster)) +
            geom_bar(position = "dodge") +
            facet_wrap(~variable,
                scales = "free", ncol = 3,
                labeller = labeller(variable = as_labeller(function(x) {
                    unlist(lapply(x, get_pretty_name))
                }))
            ) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(
                title = "Distributions of Categorical Model Variables by Cluster",
                x = NULL, y = "Count", fill = "Cluster"
            )

        print(p_cat)
        ggsave("discrepancies/2025-10-Data-Version/clustering/categorical_variable_distributions_by_cluster.pdf", plot = p_cat, width = 12, height = 8)
    }
}

# Summarized "heatmap" of variable means per cluster (numeric) and mode per cluster (categorical)
# For numeric variables: compute mean and SD per cluster, compare to overall mean/SD
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
        # Discretize z-score into categories
        cat = case_when(
            z > 1.5 ~ "++",
            z > 0.5 ~ "+",
            z > -0.5 ~ "=",
            z > -1.5 ~ "-",
            TRUE ~ "--"
        ),
        # Use SD to indicate variation (e.g., as alpha or size)
        variation = sd_val / overall_sd
    )

# For categorical variables: compute most frequent value per cluster and its proportion
cat_summary <- if (exists("df_cat")) {
    df_cat %>%
        group_by(variable, cluster, value) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(variable, cluster) %>%
        mutate(prop = n / sum(n)) %>%
        arrange(variable, cluster, desc(prop)) %>%
        slice(1) %>%
        ungroup() %>%
        mutate(
            # Use prop to indicate how dominant the mode is (variation)
            variation = 1 - prop
        )
} else {
    tibble()
}

# Special summary for TC_Path_minus_TC_AI
tc_path_summary <- NULL
if (has_tc_path) {
    tc_path_summary <- df_tc_path %>%
        group_by(variable, cluster) %>%
        summarise(
            mean_val = mean(value, na.rm = TRUE),
            sd_val = sd(value, na.rm = TRUE),
            n = n(),
            .groups = "drop"
        ) %>%
        mutate(
            # Use hard cutoffs for category
            cat = case_when(
                mean_val > 40 ~ "++",
                mean_val > 20 ~ "+",
                mean_val > 0 ~ "=",
                mean_val > -20 ~ "-",
                TRUE ~ "--"
            ),
            variation = sd_val / (max(abs(mean_val), 1e-6)) # avoid div by zero
        )
} else {
    # find the means etc for the clusters
    tc_path_summary <- df_all %>%
        group_by(cluster) %>%
        summarise(
            mean_val = mean(!!response_var, na.rm = TRUE),
            sd_val = sd(!!response_var, na.rm = TRUE),
            n = n(),
            .groups = "drop"
        ) %>%
        mutate(
            # Use hard cutoffs for category
            cat = case_when(
                mean_val > 20 ~ "++",
                mean_val > 10 ~ "+",
                mean_val < -10 ~ "-",
                mean_val < -20 ~ "--",
                TRUE ~ "="
            ),
            variation = sd_val / (max(abs(mean_val), 1e-6)) # avoid div by zero
        )
}

# Prepare numeric summary for plotting
num_plot_df <- num_summary %>%
    mutate(
        value_type = "numeric",
        display = cat, # use discretized z-score category for fill
        fill_legend = "Mean vs\nOverall",
        alpha_val = 1 - pmin(variation, 0.7)
    ) %>%
    select(cluster, variable, display, value_type, fill_legend, alpha_val)

# Prepare categorical summary for plotting
cat_plot_df <- cat_summary %>%
    mutate(
        value_type = "categorical",
        display = value, # use most frequent value for fill
        fill_legend = "Most Frequent",
        alpha_val = 1 - pmin(variation, 0.7)
    ) %>%
    select(cluster, variable, display, value_type, fill_legend, alpha_val)

# Prepare TC_Path_minus_TC_AI summary for plotting
tc_path_plot_df <- NULL
if (!is.null(tc_path_summary)) {
    tc_path_plot_df <- tc_path_summary %>%
        mutate(
            variable = "TC_Path_minus_TC_AI",
            value_type = "TC_Path_minus_TC_AI",
            display = cat,
            fill_legend = "Hard Cutoff",
            alpha_val = 1 - pmin(variation, 0.7)
        ) %>%
        select(cluster, display, mean_val, variable, value_type, fill_legend, alpha_val)

    # Ensure cluster is a factor
    if (!is.factor(tc_path_plot_df$cluster)) {
        tc_path_plot_df$cluster <- factor(tc_path_plot_df$cluster)
    }
}

# Combine all for heatmap
heatmap_df <- bind_rows(num_plot_df, cat_plot_df, tc_path_plot_df)

# Order variables: numeric first, then categorical, then TC_Path_minus_TC_AI, or keep original order
heatmap_df$variable <- factor(
    heatmap_df$variable,
    levels = unique(c(num_summary$variable, cat_summary$variable, "TC_Path_minus_TC_AI"))
)

# Plot: use fill for display, facet for value_type, and alpha for variation
library(ggplot2)
library(ggnewscale)



ggplot(heatmap_df, aes(x = factor(cluster), y = variable)) +
    # First: fill for numeric variables
    geom_tile(
        data = filter(heatmap_df, value_type == "numeric") %>% mutate(
            display = factor(display, levels = c("--", "-", "=", "+", "++"))
        ),
        aes(fill = display, alpha = alpha_val),
        color = "white", show.legend = TRUE
    ) +
    scale_fill_manual(
        name = NULL,
        values = c(
            "#b2182b",
            "#ef8a62",
            "#f7f7f7",
            "#67a9cf",
            "#2166ac"
        ),
        # breaks = c("++", "+", "=", "-", "--"),
        drop = FALSE, # ensures all breaks appear in legend even if not in data
        guide = guide_legend(order = 1, override.aes = list(alpha = 1))
    ) +
    guides(
        fill = guide_legend(
            override.aes = list(alpha = 1),
            title = NULL
        ),
        alpha = "none"
    ) +
    scale_alpha(range = c(0.75, 1), guide = "none") +
    ggnewscale::new_scale_fill() +
    # Second: fill for categorical variables (Yes/No)
    geom_tile(
        data = subset(heatmap_df, value_type == "categorical") %>% mutate(
            display = factor(display, levels = c("No", "Yes"))
        ),
        aes(fill = display, alpha = alpha_val),
        color = "white"
    ) +
    scale_fill_manual(
        name = NULL,
        values = c(
            "Yes" = "#2166ac",
            "No" = "#b2182b"
        ),
        # breaks = c("Yes", "No"),
        guide = guide_legend(order = 2),
        na.value = "grey90",
        drop = FALSE
    ) +
    ggnewscale::new_scale_fill() +
    # Third: fill for TC_Path_minus_TC_AI (numeric mean_val)
    geom_tile(
        data = subset(heatmap_df, value_type == "TC_Path_minus_TC_AI") %>% mutate(variable = "TCC Path-AI Discrepancy"),
        aes(fill = mean_val, alpha = alpha_val),
        color = "white"
    ) +
    scale_fill_gradient2(
        name = NULL,
        low = "#b2182b",
        mid = "#f7f7f7",
        high = "#2166ac",
        midpoint = 0,
        na.value = "grey90"
    ) +
    labs(
        title = "Cluster-wise summary of model variables",
        x = "Cluster",
        y = "Variable"
    ) +
    scale_y_discrete(labels = function(x) {
        sapply(x, function(lbl) {
            pretty <- get_pretty_name(lbl)
            # Insert a newline after every 7 characters, at the next space
            chars <- unlist(strsplit(pretty, ""))
            out <- ""
            count <- 0
            for (i in seq_along(chars)) {
                out <- paste0(out, chars[i])
                count <- count + 1
                if (count >= 7 && chars[i] == " ") {
                    out <- paste0(out, "\n")
                    count <- 0
                }
            }
            out
        })
    }) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 0),
        panel.spacing.y = unit(0.5, "lines")
    ) +
    coord_flip() +
    theme(legend.position = "bottom") +
    ylab(NULL)

ggsave("cluster_variable_summary_heatmap.pdf", width = 6 * 3, height = 6)