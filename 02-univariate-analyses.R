# univariate test all variables


# Helper to test one variable vs TC_Path_minus_TC_AI
univariate_test <- function(df, var, outcome = "TCC_Patho_minus_TCC_AI") {
    v <- df[[var]][!is.na(df[[var]])]
    y <- df[[outcome]][!is.na(df[[var]])]
    res <- NULL
    if (is.numeric(v)) {
        # Spearman correlation
        test <- suppressWarnings(cor.test(v, y, method = "spearman", use = "complete.obs"))
        res <- tibble(
            variable = var,
            type = "numeric",
            stat = "spearman_rho",
            value = unname(test$estimate),
            p.value = test$p.value
        )
    } else if (is.factor(v) || is.character(v) || is.ordered(v)) {
        # Pairwise Wilcoxon (vs reference)
        v_fac <- as.factor(v)
        lvls <- levels(v_fac)
        if (length(lvls) < 2) {
            return(NULL)
        }
        pairs <- combn(lvls, 2, simplify = FALSE)
        res <- map_dfr(pairs, function(pair) {
            group1 <- y[v_fac == pair[1]]
            group2 <- y[v_fac == pair[2]]
            if (length(group1) < 2 || length(group2) < 2) {
                return(NULL)
            }
            test <- suppressWarnings(wilcox.test(group1, group2))
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


# Run for all variables except the outcome itself
all_vars <- setdiff(names(data_df_pre_scaling), response_vars)



# Loop through each response variable
for (response_var in response_vars) {
    univ_results <- map_dfr(all_vars, ~ univariate_test(data_df_pre_scaling_NAd_sampletypes, .x, outcome = response_var))

    # Remove duplicate numeric variables with same value and similar names (differing only by up to 4-letter prefix)
    univ_results <- univ_results %>%
        group_by(type) %>%
        mutate(
            base_var = if_else(
                type == "numeric",
                sub("^[A-Za-z0-9_]{0,4}_", "", variable),
                variable
            )
        ) %>%
        ungroup()

    # For numeric type: keep only the row with the shortest prefix (i.e., the variable name without the prefix)
    univ_results <- univ_results %>%
        group_by(type, base_var) %>%
        filter(
            type != "numeric" |
                # For numeric: keep only the row with the variable name exactly equal to base_var
                !(type == "numeric" & n() > 1 & variable != base_var & all(value == value[1]))
        ) %>%
        ungroup() %>%
        dplyr::select(-base_var)

    univ_results$variable <- univ_results$variable %>%
        lapply(get_pretty_name) %>%
        unlist()

    univ_results <- univ_results %>%
        mutate(
            stat = case_when(
                type == "categorical" ~ {
                    parts <- strsplit(stat, " vs ")
                    out <- sapply(seq_along(parts), function(i) {
                        x <- parts[[i]]
                        if (length(x) == 2) {
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

    # Plot numeric variables: Spearman rho
    vars_to_exclude <- if (response_var == "TCC_Patho_minus_TCC_FMI") {
        c(variables_not_in_path_vs_fmi, variables_not_in_model) %>%
            lapply(get_pretty_name) %>%
            unlist()
    } else {
        variables_not_in_model %>%
            lapply(get_pretty_name) %>%
            unlist()
    }

    plot_numeric <- univ_results %>%
        filter(type == "numeric", !(variable %in% vars_to_exclude)) %>%
        mutate(variable = factor(variable,
                                 levels = variable[order(value)])) %>%
        ggplot(aes(x = value, y = variable, color = p.value < 0.05)) +
        geom_point(size = 3) +
        scale_color_manual(values = c("grey60", "firebrick")) +
        labs(
          x = "Spearman's rho", y = NULL, color = "p < 0.05",
          title = paste("Outcome:", response_var)
        ) +
        theme_bw(14)

    # Plot categorical variables: effect size (median diff)
    plot_categorical <- univ_results %>%
        filter(type == "categorical", !(variable %in% vars_to_exclude)) %>%
        mutate(variable = factor(variable)) %>%
        ggplot(aes(x = stat, y = variable, fill = p.value < 0.05)) +
        geom_tile(color = "white", height = 0.8, width = 0.8) +
        scale_fill_manual(values = c("grey90", "firebrick")) +
        labs(
            x = "Group comparison", y = NULL, fill = "p < 0.05",
            title = paste("Outcome:", response_var)
        ) +
        theme_bw(14) +
        # scale_x_discrete(labels = NULL) +
        scale_y_discrete(labels = NULL) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.ticks = element_blank(),
            # axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        facet_wrap(~variable, scales = "free", ncol = 2)

    

    dir.create(file.path(project_dir, "output/univar"), recursive = TRUE, showWarnings = FALSE)

    # Show plots with response variable in filename
    plot_numeric %>% ggsave(filename = file.path(project_dir, "output/univar", paste0("univariate_numeric_", response_var, ".pdf")), width = 10, height = 4, dpi = 300)
    plot_categorical %>% ggsave(filename = file.path(project_dir, "output/univar", paste0("univariate_categorical_", response_var, ".pdf")), width = 10, height = 30, dpi = 300)
    write_tsv(univ_results, file = file.path(project_dir, "output/univar", paste0("univariate_results_", response_var, ".tsv")))
}
