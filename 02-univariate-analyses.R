################################################################################
# Script: 02-univariate-analyses.R
# Purpose: Perform univariate statistical tests for all predictors vs TCC discrepancies
# Author: Jan-Niklas Runge
# 
# Description:
#   - Tests each predictor individually against each discrepancy outcome
#   - Numeric predictors: Spearman correlation
#   - Categorical predictors: Pairwise Wilcoxon tests between groups
#   - Generates visualization plots and exports results tables
# 
# Inputs:
#   - output/processed_data/data_df_pre_scaling.rds 
#   - output/processed_data/data_df_pre_scaling_NAd_sampletypes.rds
#   - output/processed_data/variables.rds
# 
# Outputs:
#   - output/univar/univariate_numeric_[outcome].pdf: Spearman rho plots
#   - output/univar/univariate_categorical_[outcome].pdf: Heatmaps of group comparisons
#   - output/univar/univariate_results_[outcome].tsv: Full results tables
# 

# Dependencies: See 00-universal-dependencies.R for package requirements
################################################################################

# Configuration ----------------------------------------------------------------
# Significance threshold for highlighting in plots
SIGNIFICANCE_ALPHA <- 0.05

# Plot dimensions
NUMERIC_PLOT_WIDTH <- 10
NUMERIC_PLOT_HEIGHT <- 4
CATEGORICAL_PLOT_WIDTH <- 10
CATEGORICAL_PLOT_HEIGHT <- 30
PLOT_DPI <- 300

# Minimum group size for valid statistical test
MIN_GROUP_SIZE <- 2

# Helper functions -------------------------------------------------------------

#' Perform univariate test for one predictor vs one outcome
#' 
#' @param df Data frame containing variables
#' @param var Character string, name of predictor variable
#' @param outcome Character string, name of outcome variable (default: "TCC_Patho_minus_TCC_AI")
#' @return Tibble with test results (variable, type, stat, value, p.value)
#' @details 
#'   - Numeric predictors: Spearman correlation
#'   - Categorical predictors: Pairwise Wilcoxon tests for all level combinations
#'   - Returns NULL if insufficient data for testing
univariate_test <- function(df, var, outcome = "TCC_Patho_minus_TCC_AI") {
    # Extract variable and outcome, removing NAs
    v <- df[[var]][!is.na(df[[var]])]
    y <- df[[outcome]][!is.na(df[[var]])]
    res <- NULL
    
    if (is.numeric(v)) {
        # Numeric predictor: Spearman correlation
        # Rationale: Robust to non-normality and monotonic relationships
        test <- suppressWarnings(cor.test(v, y, method = "spearman", use = "complete.obs"))
        res <- tibble(
            variable = var,
            type = "numeric",
            stat = "spearman_rho",
            value = unname(test$estimate),
            p.value = test$p.value
        )
    } else if (is.factor(v) || is.character(v) || is.ordered(v)) {
        # Categorical predictor: Pairwise Wilcoxon tests
        # Rationale: Non-parametric test for differences in median outcome between groups
        v_fac <- as.factor(v)
        lvls <- levels(v_fac)
        
        # Need at least 2 levels for comparison
        if (length(lvls) < 2) {
            return(NULL)
        }
        
        # Generate all pairwise combinations
        pairs <- combn(lvls, 2, simplify = FALSE)
        
        res <- map_dfr(pairs, function(pair) {
            group1 <- y[v_fac == pair[1]]
            group2 <- y[v_fac == pair[2]]
            
            # Require minimum group size for valid test
            if (length(group1) < MIN_GROUP_SIZE || length(group2) < MIN_GROUP_SIZE) {
                return(NULL)
            }
            
            test <- suppressWarnings(wilcox.test(group1, group2))
            
            # Effect size: difference in medians
            eff <- median(group1, na.rm = TRUE) - median(group2, na.rm = TRUE)
            
            tibble(
                variable = var,
                type = "categorical",
                stat = paste0(pair[1], " vs ", pair[2]),
                value = eff,
                p.value = test$p.value
            )
        })
    }
    res
}

# Run univariate analyses for all outcomes ------------------------------------

# Identify all predictor variables (exclude response variables)
all_vars <- setdiff(names(data_df_pre_scaling), response_vars)

# Loop through each response variable
for (response_var in response_vars) {
    message(paste("\n=== Processing outcome:", response_var, "==="))
    
    # Run univariate tests for all predictors ---------------------------------
    univ_results <- map_dfr(
        all_vars, 
        ~ univariate_test(data_df_pre_scaling_NAd_sampletypes, .x, outcome = response_var)
    )

    # Apply pretty names for presentation
    univ_results$variable <- univ_results$variable %>%
        lapply(get_pretty_name) %>%
        unlist()

    # Format categorical test labels for compact display ----------------------
    # Rationale: Shorten long factor level names to fit in plots
    univ_results <- univ_results %>%
        mutate(
            stat = case_when(
                type == "categorical" ~ {
                    parts <- strsplit(stat, " vs ")
                    out <- sapply(seq_along(parts), function(i) {
                        x <- parts[[i]]
                        if (length(x) == 2) {
                            # Abbreviate to first 3 characters
                            paste(substr(x[1], 1, 3), "vs", substr(x[2], 1, 3))
                        } else {
                            stat[i]
                        }
                    })
                    as.character(out)
                },
                TRUE ~ as.character(stat)
            )
        )

    # Determine variables to exclude from plots -------------------------------
    # Rationale: Path vs FMI comparisons exclude image quality variables
    vars_to_exclude <- if (response_var == "TCC_Patho_minus_TCC_FMI") {
        c(variables_not_in_path_vs_fmi, variables_not_in_model) %>%
            lapply(get_pretty_name) %>%
            unlist()
    } else {
        variables_not_in_model %>%
            lapply(get_pretty_name) %>%
            unlist()
    }

    # Plot numeric variables: Spearman rho ------------------------------------
    plot_numeric <- univ_results %>%
        filter(type == "numeric", !(variable %in% vars_to_exclude)) %>%
        mutate(variable = factor(variable, levels = variable[order(value)])) %>%
        ggplot(aes(x = value, y = variable, color = p.value < SIGNIFICANCE_ALPHA)) +
        geom_point(size = 3) +
        scale_color_manual(
            values = c("grey60", "firebrick"),
            labels = c("No", "Yes")
        ) +
        labs(
            x = "Spearman's rho", 
            y = NULL, 
            color = paste0("p < ", SIGNIFICANCE_ALPHA),
            title = paste("Univariate associations with", get_pretty_name(response_var))
        ) +
        theme_bw(14) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom"
        )

    # Plot categorical variables: heatmap of p-values -------------------------
    # Rationale: Heatmap allows visualization of many pairwise comparisons
    plot_categorical <- univ_results %>%
        filter(type == "categorical", !(variable %in% vars_to_exclude)) %>%
        mutate(variable = factor(variable)) %>%
        ggplot(aes(x = stat, y = variable, fill = p.value < SIGNIFICANCE_ALPHA)) +
        geom_tile(color = "white", height = 0.8, width = 0.8) +
        scale_fill_manual(
            values = c("grey90", "firebrick"),
            labels = c("No", "Yes")
        ) +
        labs(
            x = "Group comparison", 
            y = NULL, 
            fill = paste0("p < ", SIGNIFICANCE_ALPHA),
            title = paste("Pairwise comparisons for", get_pretty_name(response_var))
        ) +
        theme_bw(14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom",
            strip.background = element_rect(fill = "grey95"),
            strip.text = element_text(face = "bold")
        ) +
        facet_wrap(~variable, scales = "free", ncol = 2)

    # Save outputs -------------------------------------------------------------
    output_dir <- file.path(project_dir, "output/univar")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Save plots
    ggsave(
        filename = file.path(output_dir, paste0("univariate_numeric_", response_var, ".pdf")), 
        plot = plot_numeric,
        width = NUMERIC_PLOT_WIDTH, 
        height = NUMERIC_PLOT_HEIGHT, 
        dpi = PLOT_DPI
    )
    
    ggsave(
        filename = file.path(output_dir, paste0("univariate_categorical_", response_var, ".pdf")), 
        plot = plot_categorical,
        width = CATEGORICAL_PLOT_WIDTH, 
        height = CATEGORICAL_PLOT_HEIGHT, 
        dpi = PLOT_DPI
    )
    
    # Save results table
    write_tsv(
        univ_results, 
        file = file.path(output_dir, paste0("univariate_results_", response_var, ".tsv"))
    )
    
    message(paste("Results saved for", response_var))
}

message("\n=== Univariate analyses complete ===")
