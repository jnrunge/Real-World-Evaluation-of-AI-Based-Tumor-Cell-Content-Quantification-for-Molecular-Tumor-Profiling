# ---- Plot each variable vs TC_Path_minus_TC_AI with p-value annotation ----

source(file.path(project_dir, "00-universal-dependencies.R"))

# Loop over all response variables
for (response_var in response_vars) {
    univ_results <- read_tsv(file.path(project_dir, "output/univar", paste0("univariate_results_", response_var, ".tsv")))
    plot_vs_discrepancy <- function(df, var, univ_results, get_pretty_name) {
        if (!(get_pretty_name(var) %in% univ_results$variable)) {
            return(NULL) # Skip if variable not in univ_results
        }
        v <- df[[var]][!is.na(df[[var]])]
        y <- df[[response_var]][!is.na(df[[var]])]
        df <- df %>% filter(!is.na(.data[[var]]))
        pretty_var <- get_pretty_name(var)
        # Get p-value(s) for this variable from univ_results
        pvals <- univ_results %>% filter(variable == get_pretty_name(var))
        if (is.numeric(v)) {
            val <- pvals$value[1]
            val_label <- if (!is.na(val)) paste0("rs = ", signif(val, 2)) else ""
            # For Spearman, show a smooth trend (loess) and annotate with Spearman's rho
            rho <- pvals$value[1] %>% round(2)
            p <- ggplot(df, aes(x = .data[[var]], y = .data[[response_var]])) +
                geom_point(alpha = 0.7, color = "#1041FF", size = 2) +
                geom_smooth(method = "lm", color = "#a84632", size = 1.2, se = FALSE, linetype = "dashed") +
                labs(x = NULL, y = NULL, title = pretty_var) +
                annotate("label",
                    x = Inf, y = Inf,
                    label = paste0(val_label),
                    hjust = 1.1, vjust = 1.5, size = 5, color = "black",
                    fill = "white", label.size = 0.7
                ) +
                theme_bw(14) +
                coord_cartesian(ylim = c(-100, 100)) +
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
                geom_boxplot(fill = "#72c2ff", color = "#013d6b", alpha = 0.8, outlier.shape = NA) +
                geom_jitter(width = 0.2, alpha = 0.5, color = "#1041FF", size = 1.5) +
                labs(x = NULL, y = NULL, title = pretty_var) +
                theme_bw(14) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                coord_cartesian(ylim = c(-100, 100))
            # Add significance lines and asterisks if any
            if (!is.null(signif_df) && nrow(signif_df) > 0) {
                p <- p +
                    geom_signif(
                        comparisons = signif_df %>%
                            select(group1, group2) %>%
                            pmap(c),
                        test = "wilcox.test",
                        step_increase = 0.06,
                        vjust = 0.7,
                        textsize = 5,
                        map_signif_level = TRUE
                    )
            }
        } else {
            return(NULL)
        }
        p
    }

    # Generate plots for all variables except the outcome itself
    all_vs_discrepancy_plots <- lapply(
        setdiff(names(data_df_pre_scaling_NAd_sampletypes), response_var),
        function(var) plot_vs_discrepancy(data_df_pre_scaling_NAd_sampletypes, var, univ_results, get_pretty_name)
    )
    all_vs_discrepancy_plots <- Filter(Negate(is.null), all_vs_discrepancy_plots)

    # Combine into a grid (patchwork)
    ncol_grid <- 3
    combined_vs_discrepancy_plot <- wrap_plots(all_vs_discrepancy_plots, ncol = ncol_grid)

    # Show or save
    ggsave(file.path(project_dir, "output/univar", paste0("all_vs_", response_var, ".pdf")), combined_vs_discrepancy_plot, width = 16, height = 40, dpi = 300)



    ## manually selected plots by Mariam
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
        function(var) plot_vs_discrepancy(data_df_pre_scaling_NAd_sampletypes, var, univ_results, get_pretty_name)
    )
    manual_discrepancy_plots_list <- Filter(Negate(is.null), manual_discrepancy_plots_list)
    n_plots <- length(manual_discrepancy_plots_list)
    ncol_grid <- 2

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

    ggsave(file.path(project_dir, "output/univar", paste0("manual_discrepancy_plots_", response_var, ".pdf")), combined_manual_discrepancy_plots, width = 10, height = 18, dpi = 300)


    ai_vars_to_show <- c(
        # "Foundation.Model...TME.v1.0.1.alpha...Key.Result...Total.area.of.tissue...mm.2", -- no longer wanted; intended?
        "AI_total_stroma_area_%",
        "AI_total_cell_density_mm2",
        "AI_total_necrosis_area_%"
    )

    ai_vars_plots_list <- lapply(
        ai_vars_to_show,
        function(var) plot_vs_discrepancy(data_df_pre_scaling_NAd_sampletypes, var, univ_results, get_pretty_name)
    )
    ai_vars_plots_list <- Filter(Negate(is.null), ai_vars_plots_list)
    n_plots <- length(ai_vars_plots_list)
    ncol_grid <- 2

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

    ggsave(file.path(project_dir, "output/univar", paste0("ai_vars_vs_", response_var, ".pdf")), combined_ai_vars_plots, width = 10, height = 9, dpi = 300)

    # --- Only variables with at least one significant p-value ---
    sig_vars <- univ_results %>%
        group_by(variable) %>%
        filter(any(p.value < 0.05, na.rm = TRUE)) %>%
        pull(variable) %>%
        unique()

    # Helper to classify pretty variable names as Path or AI
    is_ai_var <- function(pretty_var) {
        grepl("^[(]?AI", pretty_var)
    }
    is_path_var <- function(pretty_var) {
        !grepl("^[(]?AI", pretty_var)
    }

    # Path variables: pretty name does NOT start with (AI)
    sig_vs_discrepancy_plots_path <- lapply(
        setdiff(names(data_df_pre_scaling_NAd_sampletypes), response_vars),
        function(var) {
            pretty_var <- get_pretty_name(var)
            if (pretty_var %in% sig_vars && is_path_var(pretty_var)) {
                plot_vs_discrepancy(data_df_pre_scaling_NAd_sampletypes, var, univ_results, get_pretty_name)
            } else {
                NULL
            }
        }
    )
    sig_vs_discrepancy_plots_path <- Filter(Negate(is.null), sig_vs_discrepancy_plots_path)

    combined_sig_vs_discrepancy_plot_path <- wrap_plots(sig_vs_discrepancy_plots_path, ncol = ncol_grid)

    ggsave(file.path(project_dir, "output/univar", paste0("sig_vs_path_", response_var, ".pdf")), combined_sig_vs_discrepancy_plot_path, width = 16, height = 20, dpi = 300)

    # AI variables: pretty name starts with AI
    sig_vs_discrepancy_plots_ai <- lapply(
        setdiff(names(data_df_pre_scaling_NAd_sampletypes), response_vars),
        function(var) {
            pretty_var <- get_pretty_name(var)
            print(pretty_var)
            if (pretty_var %in% sig_vars && is_ai_var(pretty_var)) {
                plot_vs_discrepancy(data_df_pre_scaling_NAd_sampletypes, var, univ_results, get_pretty_name)
            } else {
                NULL
            }
        }
    )
    sig_vs_discrepancy_plots_ai <- Filter(Negate(is.null), sig_vs_discrepancy_plots_ai)

    combined_sig_vs_discrepancy_plot_ai <- wrap_plots(sig_vs_discrepancy_plots_ai, ncol = ncol_grid)

    ggsave(file.path(project_dir, "output/univar", paste0("sig_vs_ai_", response_var, ".pdf")), combined_sig_vs_discrepancy_plot_ai, width = 16, height = 10, dpi = 300)

} # End of loop over response_vars
