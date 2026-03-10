################################################################################
# Script: 02c-other-univ-plots.R
# Purpose: Generate scatter/boxplot visualizations for univariate associations
# Author: Jan-Niklas Runge
# 
# Description:
#   - Creates individual plots for each predictor vs each discrepancy outcome
#   - Numeric predictors: scatter plots with linear regression and Spearman rho
#   - Categorical predictors: boxplots with significance markers for pairwise tests
# 
# 
# Outputs:
#   - output/univar/all_vs_[outcome].pdf: All predictors in grid layout
#   - output/univar/manual_discrepancy_plots_[outcome].pdf: Manually selected variables
#   - output/univar/ai_vars_vs_[outcome].pdf: AI-specific variables
#   - output/univar/sig_vs_path_[outcome].pdf: Significant pathology variables
#   - output/univar/sig_vs_ai_[outcome].pdf: Significant AI variables
# 
# Dependencies: See 00-universal-dependencies.R for package requirements
################################################################################

# Configuration ----------------------------------------------------------------
# Grid layout parameters
GRID_NCOL_ALL <- 3         # Columns for comprehensive grid
GRID_NCOL_SUBSET <- 2      # Columns for filtered subsets

# Plot dimensions (inches)
PLOT_WIDTH_ALL <- 16
PLOT_HEIGHT_ALL <- 40
PLOT_WIDTH_MANUAL <- 10
PLOT_HEIGHT_MANUAL <- 18
PLOT_WIDTH_AI <- 10
PLOT_HEIGHT_AI <- 9
PLOT_WIDTH_SIG_PATH <- 16
PLOT_HEIGHT_SIG_PATH <- 20
PLOT_WIDTH_SIG_AI <- 16
PLOT_HEIGHT_SIG_AI <- 10
PLOT_DPI <- 300

# Y-axis limits for consistency across plots
YLIM_DISCREPANCY <- c(-100, 100)

# Annotation positioning (relative to plot area)
ANNOT_HJUST <- 1.5
ANNOT_VJUST <-2
ANNOT_SIZE <- 5

# Point styling
POINT_ALPHA <- 0.7
POINT_COLOR <- "#1041FF"
POINT_SIZE <- 0.5

# Correlation label styling
CORR_LABEL_TEXT_SIZE <- 2
CORR_LABEL_BORDER_SIZE <- 0.25

# Line styling
LINE_COLOR <- "#a84632"
LINE_SIZE <- 0.5

# Boxplot styling
BOX_FILL <- "#72c2ff"
BOX_COLOR <- "#013d6b"
BOX_ALPHA <- 0.8
BOX_LINE_SIZE <- 0.3
JITTER_WIDTH <- 0.2
JITTER_ALPHA <- 0.5
JITTER_SIZE <- 0.5

# Significance annotation
SIGNIF_STEP_INCREASE <- 0.08
SIGNIF_VJUST <- 0.7
SIGNIF_TEXTSIZE <- 3

# Load dependencies ------------------------------------------------------------
source(file.path(project_dir, "00-universal-dependencies.R"))

# Helper functions -------------------------------------------------------------
variable_glossary <- pretty_names
#' Get short display name for variable
#' @param var_name Character string, variable name
#' @return Character string (short name) or NA if not available
get_short_name <- function(var_name) {
    if (exists("variable_glossary") && "pretty_name_short" %in% names(variable_glossary)) {
        val <- variable_glossary$pretty_name_short[match(var_name, variable_glossary$variable)]
        if (length(val) > 0) return(val)
    }
    return(NA_character_)
}

#' Create plot for a single predictor vs outcome
#' 
#' @param df Data frame containing variables
#' @param var Character string, name of predictor variable
#' @param univ_results Tibble with univariate test results
#' @param get_pretty_name Function to convert variable names to display names
#' @param get_short_name Function to convert variable names to short display names (optional)
#' @return ggplot object or NULL if variable not in results
#' @details 
#'   - Numeric: scatter plot with linear trend and Spearman rho annotation
#'   - Categorical: boxplot with significance brackets for p < 0.05 comparisons
plot_vs_discrepancy <- function(df, var, univ_results, get_pretty_name, get_short_name = NULL) {
    if (!(get_pretty_name(var) %in% univ_results$variable)) {
        return(NULL) # Skip if variable not in univ_results
    }
    v <- df[[var]][!is.na(df[[var]])]
    y <- df[[response_var]][!is.na(df[[var]])]
    df <- df %>% filter(!is.na(.data[[var]]))
    pretty_var <- get_pretty_name(var)
    title_text_size <- 6

    # Determine plot title: prefer short name if available and not empty
    plot_title <- pretty_var
    short_name <- NA
    
    if (!is.null(get_short_name)) {
        short_name <- get_short_name(var)
    }
    
    # Fallback: if short name lookup failed by variable key (e.g. due to dots vs spaces),
    # try looking up using the pretty_name which we already resolved successfully
    if ((length(short_name) == 0 || is.na(short_name) || short_name == "") && 
        exists("variable_glossary") && "pretty_name" %in% names(variable_glossary)) {
        idx <- match(pretty_var, variable_glossary$pretty_name)
        if (!is.na(idx)) {
            val <- variable_glossary$pretty_name_short[idx]
            if (!is.na(val) && val != "") short_name <- val
        }
    }

    if (length(short_name) > 0 && !is.na(short_name) && short_name != "") {
        plot_title <- short_name
    }

    # Get p-value(s) for this variable from univ_results
    pvals <- univ_results %>% filter(variable == get_pretty_name(var))
    if (is.numeric(v)) {
        val <- pvals$value[1]
        val_label <- if (!is.na(val)) paste0("rs = ", signif(val, 2)) else ""
        # For Spearman, show a smooth trend (loess) and annotate with Spearman's rho
        rho <- pvals$value[1] %>% round(2)
        p <- ggplot(df, aes(x = .data[[var]], y = .data[[response_var]])) +
            geom_point(alpha = POINT_ALPHA, color = POINT_COLOR, size = POINT_SIZE) +
            geom_smooth(method = "lm", color = LINE_COLOR, size = LINE_SIZE, se = FALSE, linetype = "dashed") +
            labs(x = NULL, y = NULL, title = plot_title) +
            annotate("label",
                x = Inf, y = Inf,
                label = paste0(val_label),
                hjust = ANNOT_HJUST, vjust = ANNOT_VJUST,
                size = CORR_LABEL_TEXT_SIZE, color = "black",
                fill = "white", label.size = CORR_LABEL_BORDER_SIZE
            ) +
            theme_bw(8) +
            theme(plot.title = element_text(size = title_text_size)) +
            coord_cartesian(ylim = YLIM_DISCREPANCY) +
            scale_x_continuous(labels = scales::comma)
    } else if (is.factor(v) || is.character(v) || is.ordered(v)) {
        v_fac <- if (is.factor(v)) v else factor(v)
        # For categorical, show boxplot with significance asterisks for significant pairwise comparisons
        # Prepare pairwise comparisons and significance
        pvals_sig <- pvals %>% filter(p.value < 0.05)
        # Extract pairs from stat (format: "ABC vs DEF")
        if (nrow(pvals_sig) > 0) {
            comparisons <- strsplit(pvals_sig$stat, " vs ")
            # Get the actual factor levels (first three letters matched to levels)
            lvl_map <- setNames(levels(v_fac), substr(levels(v_fac), 1, 3))
            comp_list <- lapply(comparisons, function(pair) {
                # Find full level names by matching first 3 letters
                l1 <- lvl_map[[pair[1]]]
                l2 <- lvl_map[[pair[2]]]
                c(l1, l2)
            })
            # Prepare annotation data frame for ggsignif
            signif_df <- tibble::tibble(
                group1 = sapply(comp_list, `[`, 1),
                group2 = sapply(comp_list, `[`, 2),
                y_position = max(y, na.rm = TRUE) + seq(5, 20, length.out = length(comp_list)),
                label = "*"
            )
        } else {
            signif_df <- NULL
        }
        p <- ggplot(df, aes(x = v_fac, y = .data[[response_var]])) +
            geom_boxplot(
                fill = BOX_FILL,
                color = BOX_COLOR,
                alpha = BOX_ALPHA,
                outlier.shape = NA,
                size = BOX_LINE_SIZE
            ) +
            geom_jitter(width = JITTER_WIDTH, alpha = JITTER_ALPHA, color = POINT_COLOR, size = JITTER_SIZE) +
            labs(x = NULL, y = NULL, title = plot_title) +
            theme_bw(8) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(size = title_text_size)
            ) +
            coord_cartesian(ylim = YLIM_DISCREPANCY)
        # Add significance lines and asterisks if any
        if (!is.null(signif_df) && nrow(signif_df) > 0) {
            p <- p +
                geom_signif(
                  comparisons = signif_df %>%
                    select(group1, group2) %>%
                    pmap(c),
                  test = "wilcox.test",
                  step_increase = SIGNIF_STEP_INCREASE,
                  vjust = SIGNIF_VJUST,
                  textsize = SIGNIF_TEXTSIZE,
                  map_signif_level = TRUE,
                  size = 0.3
                )
        }
    } else {
        return(NULL)
    }
    p
}

#' Classify variable as AI-sourced based on pretty name
#' @param pretty_var Character string, formatted variable name
#' @return Logical, TRUE if AI variable
is_ai_var <- function(pretty_var) {
    grepl("^[(]?AI", pretty_var)
}

#' Classify variable as pathology-sourced (non-AI)
#' @param pretty_var Character string, formatted variable name
#' @return Logical, TRUE if pathology variable
is_path_var <- function(pretty_var) {
    !grepl("^[(]?AI", pretty_var)
}

# Main analysis loop -----------------------------------------------------------
# Loop over all response variables
for (response_var in response_vars) {
    message(paste("\n=== Generating plots for outcome:", response_var, "==="))
    
    # Load univariate results --------------------------------------------------
    univ_results <- read_tsv(file.path(project_dir, "output/univar", paste0("univariate_results_", response_var, ".tsv")))
    
    # Generate plots for all variables -----------------------------------------
    # ---- Plot each variable vs TC_Path_minus_TC_AI with p-value annotation ----

    # Generate plots for all variables except the outcome itself
    all_vs_discrepancy_plots <- lapply(
        setdiff(names(data_df_pre_scaling_NAd_sampletypes), response_var),
        function(var) plot_vs_discrepancy(data_df_pre_scaling_NAd_sampletypes, var, univ_results, get_pretty_name, get_short_name)
    )
    all_vs_discrepancy_plots <- Filter(Negate(is.null), all_vs_discrepancy_plots)

    # Combine into a grid (patchwork) and save in chunks of 10
    ncol_grid <- GRID_NCOL_ALL
    plots_per_file <- 9
    plot_chunks <- split(all_vs_discrepancy_plots, ceiling(seq_along(all_vs_discrepancy_plots) / plots_per_file))

    for (i in seq_along(plot_chunks)) {
        chunk_plots <- plot_chunks[[i]]
        combined_vs_discrepancy_plot <- wrap_plots(chunk_plots, ncol = ncol_grid)
        
        # Calculate height dynamically based on number of rows (approx 5 inches per row)
        n_rows <- ceiling(length(chunk_plots) / ncol_grid)
        chunk_height <- n_rows * 4

        ggsave(file.path(project_dir, "output/univar", paste0("all_vs_", response_var, "_part", i, ".pdf")), 
               combined_vs_discrepancy_plot, width = 12, height = chunk_height, units="cm", dpi = PLOT_DPI)
    }

    # Manually selected plots --------------------------------------------------
    # Rationale: Key variables identified by clinical experts (Mariam)
    manual_vars_to_show <- c(
        names(data_df_pre_scaling_NAd_sampletypes)[grepl("Microto", names(data_df_pre_scaling_NAd_sampletypes))],
        names(data_df_pre_scaling_NAd_sampletypes)[grepl("Sample.*typ", names(data_df_pre_scaling_NAd_sampletypes))],
        names(data_df_pre_scaling_NAd_sampletypes)[grepl("Immu.*agg", names(data_df_pre_scaling_NAd_sampletypes))],
        names(data_df_pre_scaling_NAd_sampletypes)[grepl("Necro", names(data_df_pre_scaling_NAd_sampletypes))],
        names(data_df_pre_scaling_NAd_sampletypes)[grepl("ructur", names(data_df_pre_scaling_NAd_sampletypes))],
        names(data_df_pre_scaling_NAd_sampletypes)[grepl("Acute", names(data_df_pre_scaling_NAd_sampletypes))],
        names(data_df_pre_scaling_NAd_sampletypes)[grepl("Path_evaluable_cancer_area_marked_as_artifact_by_AI", names(data_df_pre_scaling_NAd_sampletypes))]
    )

    manual_discrepancy_plots_list <- lapply(
        manual_vars_to_show,
        function(var) plot_vs_discrepancy(data_df_pre_scaling_NAd_sampletypes, var, univ_results, get_pretty_name, get_short_name)
    )
    manual_discrepancy_plots_list <- Filter(Negate(is.null), manual_discrepancy_plots_list)
    n_plots <- length(manual_discrepancy_plots_list)
    ncol_grid <- GRID_NCOL_SUBSET

    # Handle odd number of plots: last plot spans full width
    if (n_plots %% ncol_grid == 1 && n_plots > 1) {
        # Odd number: wrap all but last, last fills full width
        top_plots <- manual_discrepancy_plots_list[1:(n_plots - 1)]
        last_plot <- manual_discrepancy_plots_list[[n_plots]]
        combined_manual_discrepancy_plots <- wrap_plots(top_plots, ncol = ncol_grid) / last_plot +
            plot_layout(heights = c((n_plots - 1) / ncol_grid, 1))
    } else {
        # Even number: wrap all
        combined_manual_discrepancy_plots <- wrap_plots(manual_discrepancy_plots_list, ncol = ncol_grid)
    }

    ggsave(file.path(project_dir, "output/univar", paste0("manual_discrepancy_plots_", response_var, ".pdf")), combined_manual_discrepancy_plots, width = PLOT_WIDTH_MANUAL, height = PLOT_HEIGHT_MANUAL, dpi = PLOT_DPI)

    # AI-specific variables ----------------------------------------------------
    # Rationale: Core AI-derived features for direct AI vs Path comparison
    ai_vars_to_show <- c(
        # "Foundation.Model...TME.v1.0.1.alpha...Key.Result...Total.area.of.tissue...mm.2", -- no longer wanted; intended?
        "AI_total_stroma_area_%",
        "AI_total_cell_density_mm2",
        "AI_total_necrosis_area_%"
    )

    ai_vars_plots_list <- lapply(
        ai_vars_to_show,
        function(var) plot_vs_discrepancy(data_df_pre_scaling_NAd_sampletypes, var, univ_results, get_pretty_name, get_short_name)
    )
    ai_vars_plots_list <- Filter(Negate(is.null), ai_vars_plots_list)
    n_plots <- length(ai_vars_plots_list)
    ncol_grid <- GRID_NCOL_SUBSET

    if (n_plots %% ncol_grid == 1 && n_plots > 1) {
        # Odd number: wrap all but last, last fills full width
        top_plots <- ai_vars_plots_list[1:(n_plots - 1)]
        last_plot <- ai_vars_plots_list[[n_plots]]
        combined_ai_vars_plots <- wrap_plots(top_plots, ncol = ncol_grid) / last_plot +
            plot_layout(heights = c((n_plots - 1) / ncol_grid, 1))
    } else {
        # Even number: wrap all
        combined_ai_vars_plots <- wrap_plots(ai_vars_plots_list, ncol = ncol_grid)
    }

    ggsave(file.path(project_dir, "output/univar", paste0("ai_vars_vs_", response_var, ".pdf")), combined_ai_vars_plots, width = PLOT_WIDTH_AI, height = PLOT_HEIGHT_AI, dpi = PLOT_DPI)

    # Significant variables only -----------------------------------------------
    # Rationale: Focus on statistically significant associations (p < 0.05)
    sig_vars <- univ_results %>%
        group_by(variable) %>%
        filter(any(p.value < 0.05, na.rm = TRUE)) %>%
        pull(variable) %>%
        unique()

    # Path variables: pretty name does NOT start with (AI)
    sig_vs_discrepancy_plots_path <- lapply(
        setdiff(names(data_df_pre_scaling_NAd_sampletypes), response_vars),
        function(var) {
            pretty_var <- get_pretty_name(var)
            if (pretty_var %in% sig_vars && is_path_var(pretty_var)) {
                plot_vs_discrepancy(data_df_pre_scaling_NAd_sampletypes, var, univ_results, get_pretty_name, get_short_name)
            } else {
                NULL
            }
        }
    )
    sig_vs_discrepancy_plots_path <- Filter(Negate(is.null), sig_vs_discrepancy_plots_path)

    combined_sig_vs_discrepancy_plot_path <- wrap_plots(sig_vs_discrepancy_plots_path, ncol = ncol_grid)

    ggsave(file.path(project_dir, "output/univar", paste0("sig_vs_path_", response_var, ".pdf")), combined_sig_vs_discrepancy_plot_path, width = PLOT_WIDTH_SIG_PATH, height = PLOT_HEIGHT_SIG_PATH, dpi = PLOT_DPI)

    # AI variables: pretty name starts with AI
    sig_vs_discrepancy_plots_ai <- lapply(
        setdiff(names(data_df_pre_scaling_NAd_sampletypes), response_vars),
        function(var) {
            pretty_var <- get_pretty_name(var)
            print(pretty_var)
            if (pretty_var %in% sig_vars && is_ai_var(pretty_var)) {
                plot_vs_discrepancy(data_df_pre_scaling_NAd_sampletypes, var, univ_results, get_pretty_name, get_short_name)
            } else {
                NULL
            }
        }
    )
    sig_vs_discrepancy_plots_ai <- Filter(Negate(is.null), sig_vs_discrepancy_plots_ai)

    combined_sig_vs_discrepancy_plot_ai <- wrap_plots(sig_vs_discrepancy_plots_ai, ncol = ncol_grid)

    ggsave(file.path(project_dir, "output/univar", paste0("sig_vs_ai_", response_var, ".pdf")), combined_sig_vs_discrepancy_plot_ai, width = PLOT_WIDTH_SIG_AI, height = PLOT_HEIGHT_SIG_AI, dpi = PLOT_DPI)

    message(paste("Plots saved for", response_var))
} # End of loop over response_vars

message("\n=== Univariate plot generation complete ===")
