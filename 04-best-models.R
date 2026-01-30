################################################################################
# Script: 04-best-models.R
# Purpose: Build and evaluate best linear models for TCC discrepancy prediction
# Author: Jan-Niklas Runge
# 
# Description:
#   - Performs forward stepwise selection with multiple random initializations
#   - Selects best model based on AIC across 32 independent runs
#   - Generates predictor presence/absence heatmaps to visualize model stability
# 
# 
# Outputs:
#   - output/processed_data/best_models_classic_no_forced_interactions.rds: Best models
#   - output/model_plots/predictor_presence_[outcome].pdf: Predictor heatmaps with dendrograms
#   - output/model_plots/aic_progression_[outcome].pdf: AIC evolution plots
# 
################################################################################

# Configuration ----------------------------------------------------------------
# Number of independent stepwise runs per outcome
# Rationale: Multiple initializations reduce risk of local optima
N_REPS <- 32

# Maximum stepwise iterations per run
# Rationale: Ensures convergence of forward selection
MAX_STEPS <- 5000

# Random seed for reproducibility
RANDOM_SEED <- 1

# Plot dimensions (inches)
PLOT_WIDTH <- 12
PLOT_HEIGHT <- 10
PLOT_DPI <- 300

# Dendrogram height ratio in combined plot
DEND_HEIGHT_RATIO <- 1
MAIN_HEIGHT_RATIO <- 4

# Load dependencies ------------------------------------------------------------
pretty_names <- read_csv(file.path(project_dir, "input/variable_pretty_names.csv"))
source(file.path(project_dir, "functions/best_models.R"))

# Extract response variables ---------------------------------------------------
response_vars <- variables %>%
    filter(type == "response") %>%
    pull(variable)

# Apply variable exclusions ----------------------------------------------------
# Rationale: Remove variables excluded based on univariate analyses
variables <- variables %>%
    filter(!(variable %in% variables_not_in_model))

# Load or build best models ----------------------------------------------------
rds_file <- file.path(project_dir, "output/processed_data/best_models_classic_no_forced_interactions.rds")

if (file.exists(rds_file)) {
    message("Loading existing best models from RDS file...")
    classic_no_forced_interactions <- readRDS(rds_file)
} else {
    message("No existing models found. Running forward selection for all outcomes...")
    classic_no_forced_interactions <- list()
    
    # Run forward selection for each response variable ------------------------
    for (response_var in response_vars) {
        message(paste("\n=== Building models for outcome:", response_var, "==="))
        
        # Create variable table for this response
        var_table <- variables %>%
            filter(type != "response" | variable == response_var)
        
        # Apply Path vs FMI exclusions if applicable
        # Rationale: Image quality metrics irrelevant for FMI comparisons
        if (response_var == "TCC_Patho_minus_TCC_FMI") {
            var_table <- var_table %>%
                filter(!(variable %in% variables_not_in_path_vs_fmi))
            message("Applied Path vs FMI exclusions")
        }
        
        # Perform forward selection with multiple initializations
        classic_no_forced_interactions[[response_var]] <- perform_forward_selection(
            var_table %>% filter(type == "response" | !grepl("TCC", variable)), 
            data_df, 
            force_interaction = FALSE, 
            steps = MAX_STEPS, 
            seed = RANDOM_SEED, 
            reps = N_REPS
        )
        
        message(paste("Best model AIC:", 
                     min(sapply(classic_no_forced_interactions[[response_var]]$fits, 
                               function(x) AIC(x$model)))))
    }

    # Create output directories ------------------------------------------------
    dir.create(file.path(project_dir, "output/model_plots"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(project_dir, "output/processed_data"), showWarnings = FALSE, recursive = TRUE)
    
    # Save best models
    saveRDS(classic_no_forced_interactions, file = rds_file)
    message("Best models saved to:", rds_file)

    # Generate model comparison visualizations --------------------------------
    message("\n=== Generating model comparison plots ===")
    
    predictor_presence <- list()

    for (response_var in response_vars) {
        message(paste("Processing outcome:", response_var))
        
        # Extract all predictors from the 32 fits -----------------------------
        all_predictors <- lapply(1:N_REPS, function(i) {
            model <- classic_no_forced_interactions[[response_var]]$fits[[i]]$model
            # Get predictor names (excluding intercept and response variable)
            pred_names <- names(coef(model))[-1]
            pred_names
        })

        # Identify best fit (lowest AIC)
        chosen_fit <- classic_no_forced_interactions[[response_var]]$best_idx

        # Get unique predictors across all fits
        unique_predictors <- unique(unlist(all_predictors))

        # Create presence/absence matrix ---------------------------------------
        # Rationale: Binary matrix showing which predictors appear in which models
        presence_matrix <- matrix(0, nrow = length(unique_predictors), ncol = N_REPS)
        rownames(presence_matrix) <- unique_predictors
        colnames(presence_matrix) <- paste0("Fit_", 1:N_REPS)

        for (i in 1:N_REPS) {
            presence_matrix[all_predictors[[i]], i] <- 1
        }

        predictor_presence[[response_var]] <- as.data.frame(presence_matrix)

        # Perform hierarchical clustering on fits -----------------------------
        # Rationale: Group models with similar predictor sets
        # Calculate distance between fits based on predictor presence
        fit_dist <- dist(t(presence_matrix), method = "binary")
        fit_hclust <- hclust(fit_dist, method = "complete")

        # Get dendrogram order
        fit_order <- fit_hclust$order

        # Identify groups of identical fits -----------------------------------
        # Rationale: Multiple initializations may converge to identical models
        fit_groups <- character(N_REPS)
        group_id <- 1
        for (i in 1:N_REPS) {
            if (fit_groups[i] == "") {
                # This fit hasn't been assigned to a group yet
                fit_groups[i] <- as.character(group_id)
                # Find all other fits with identical predictor sets
                for (j in (i + 1):N_REPS) {
                    if (j <= N_REPS && all(presence_matrix[, i] == presence_matrix[, j])) {
                        fit_groups[j] <- as.character(group_id)
                    }
                }
                group_id <- group_id + 1
            }
        }

        # Create data frame for plotting with group information ---------------
        plot_data <- predictor_presence[[response_var]] %>%
            mutate(Predictor = rownames(.)) %>%
            pivot_longer(cols = starts_with("Fit_"), names_to = "Fit", values_to = "Present") %>%
            mutate(
                Fit_num = as.numeric(gsub("Fit_", "", Fit)),
                Fit = factor(Fit, levels = paste0("Fit_", fit_order)),
                is_chosen = (Fit == paste0("Fit_", chosen_fit)),
                fit_group = fit_groups[Fit_num]
            )

        # Create dendrogram plot -----------------------------------------------
        dend_data <- dendro_data(fit_hclust)
        dend_plot <- ggplot(segment(dend_data)) +
            geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
            theme_void() +
            theme(plot.margin = margin(5, 5, 0, 5))

        # Convert predictor names to pretty display names ---------------------
        # Rationale: Handles interaction terms and applies name mappings
        plot_data$Predictor <- sapply(plot_data$Predictor, function(var) {
            parts <- strsplit(var, ":", fixed = TRUE)[[1]]
            cleaned <- sapply(parts, function(p) {
                # Fuzzy match to handle minor name variations
                hits <- agrep(p, names(data_df), value = TRUE, max.distance = 20)
                if (length(hits) > 0) {
                    d <- adist(p, hits)
                    best <- hits[which.min(d)][1]
                } else {
                    best <- p
                }
                clean_var_names(best)
            })
            paste(sort(cleaned), collapse = ":")
        }, USE.NAMES = FALSE) %>%
            lapply(get_pretty_name) %>%
            unlist()

        # Plot presence/absence heatmap ----------------------------------------
        # Rationale: Visualize predictor stability across model runs
        # - White tiles: predictor absent
        # - Colored tiles: predictor present (color indicates model group)
        # - Red border: best model (lowest AIC)
        main_plot <- plot_data %>%
            ggplot(aes(x = Fit, y = Predictor, fill = interaction(factor(Present), fit_group))) +
            geom_tile(aes(color = is_chosen, linewidth = is_chosen)) +
            scale_fill_manual(values = c(
                setNames(
                    rep("white", length(unique(fit_groups))),
                    paste0("0.", unique(fit_groups))
                ),
                setNames(
                    rainbow(length(unique(fit_groups)), s = 0.7, v = 0.8),
                    paste0("1.", unique(fit_groups))
                )
            ), guide = "none") +
            scale_color_manual(values = c("FALSE" = "white", "TRUE" = "red"), guide = "none") +
            scale_linewidth_manual(values = c("FALSE" = 0.5, "TRUE" = 2), guide = "none") +
            labs(
                title = paste0("Predictor Presence Across Model Fits:\n", response_var),
                subtitle = paste0(
                    "Red border indicates chosen fit: ", chosen_fit,
                    "\nIdentical fits have same color\n",
                    length(unique(fit_groups)), " unique model(s)"
                ),
                x = "Model Fit (ordered by similarity)", y = "Predictor", fill = ""
            ) +
            theme_minimal(14) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

        # Combine dendrogram and main plot -------------------------------------
        combined_plot <- dend_plot / main_plot +
            plot_layout(heights = c(DEND_HEIGHT_RATIO, MAIN_HEIGHT_RATIO))

        # Save combined plot
        ggsave(
            file.path(project_dir, "output/model_plots", paste0("predictor_presence_", response_var, ".pdf")),
            plot = combined_plot,
            width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI
        )
        
        message(paste("Saved predictor presence plot for", response_var))
    }
}

# Generate AIC progression plots -----------------------------------------------
message("\n=== Generating AIC progression plots ===")

for (response_var in names(classic_no_forced_interactions)) {
    fits <- classic_no_forced_interactions[[response_var]]$fits

    # Combine stepwise progression data from all fits -------------------------
    # Rationale: Visualize how AIC evolves during forward selection
    steps_combined <- bind_rows(
        lapply(seq_along(fits), function(i) {
            steps_tbl <- tryCatch(
                {
                    # Handle different possible list structures
                    if (!is.null(names(fits[[i]])) && "steps" %in% names(fits[[i]])) {
                        fits[[i]]$steps
                    } else {
                        fits[[i]][[2]]
                    }
                },
                error = function(e) NULL
            )
            if (is.null(steps_tbl)) return(NULL)
            steps_tbl %>% mutate(fit_index = i)
        })
    )

    if (is.null(steps_combined) || nrow(steps_combined) == 0) {
        message(paste("No AIC progression data available for", response_var))
        next
    }

    # Create faceted plot showing AIC evolution for each fit ------------------
    p_steps <- ggplot(steps_combined, aes(x = step, y = AIC)) +
        geom_line() +
        geom_point(size = 0.5) +
        facet_wrap(~ fit_index, scales = "free_y") +
        labs(
            title = paste("AIC progression across stepwise selection:", response_var),
            x = "Step",
            y = "AIC"
        ) +
        theme_bw(14)

    ggsave(
        file.path(project_dir, "output/model_plots", paste0("aic_progression_", response_var, ".pdf")),
        plot = p_steps,
        width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI
    )
    
    message(paste("Saved AIC progression plot for", response_var))
}

message("\n=== Best model selection complete ===")
message("Models saved to:", rds_file)
message("Plots saved to:", file.path(project_dir, "output/model_plots"))

# Future work: Term importance analysis ----------------------------------------
# Rationale: Use compare_model_drop_terms(model, terms) to assess 
# significance and contribution of individual predictors in final models
# TODO: Implement systematic term importance evaluation