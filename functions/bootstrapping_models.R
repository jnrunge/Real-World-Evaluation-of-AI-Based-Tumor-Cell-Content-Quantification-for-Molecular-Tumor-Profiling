bootstrap_models <- function(data_df,variables, n_boot, response_var, seed=1337, steps=5000) {
    foreach(
        i = 1:n_boot,
        .packages = c("tidyverse", "stats")
    ) %dopar% {
        source(file.path(project_dir, "functions/manual-stepwise.R"))
        # sample by indices so we can keep them
        set.seed(seed+i)
        repeat {
            boot_idx <- sample.int(nrow(data_df), size = nrow(data_df), replace = TRUE)
            boot_sample <- data_df[boot_idx, , drop = FALSE]
            # skip degenerate bootstrap where any column is constant
            if (!any(sapply(boot_sample, function(x) length(unique(x)) == 1))) break
        }

        # only send one reponse variable
        var_table <- variables %>%
            filter(type != "response" | variable == response_var)

        boot_forward <- stepwise_aic(boot_sample, var_table, "both", trace=FALSE,enforce_interactions=FALSE, steps = steps, seed = seed + i)

        r2 <- summary(boot_forward[[1]])$r.squared
        if (r2 == 1) {
            return(NULL)
        }

        # augment return with bootstrap indices (keep duplicates) and unique rows
        # preserve original structure so boot_forward[[1]] is still the model
        boot_forward_aug <- c(
            boot_forward,
            list(
                bootstrap_idx = boot_idx,
                bootstrap_rows_unique = sort(unique(boot_idx))
            )
        )
        # also attach as attributes on the model object
        attr(boot_forward_aug[[1]], "bootstrap_idx") <- boot_idx
        attr(boot_forward_aug[[1]], "bootstrap_rows_unique") <- sort(unique(boot_idx))

        boot_forward_aug
    }
}




# Function to process bootstrap results and create selection frequency table
process_bootstrap_results <- function(boot_forward_list, data_df, n_boot) {
    selected_counts <- lapply(boot_forward_list, function(boot_forward) {
        if (is.null(boot_forward[[1]])) {
            return(NULL)
        }
        names(coef(boot_forward[[1]]))[-1]
    })
    
    n_boot_real <- sum(!unlist(lapply(selected_counts, identical, y = NULL)))
    all_sel <- unlist(selected_counts)
    selection_freq <- as.data.frame(table(all_sel)) %>%
        mutate(all_sel = as.character(all_sel)) %>%
        arrange(desc(Freq)) %>%
        filter(!(all_sel %in% c("fitted", "residuals")))
    
    # Normalize interaction terms
    selection_freq <- selection_freq %>%
        mutate(
            key = ifelse(
                str_detect(all_sel, ":"),
                sapply(strsplit(all_sel, ":"), function(x) paste(sort(x), collapse = ":")),
                all_sel
            )
        ) %>%
        group_by(key) %>%
        summarise(Freq = sum(Freq), .groups = "drop") %>%
        rename(all_sel = key) %>%
        arrange(desc(Freq))
    
    # Clean variable names
    selection_freq$all_sel <- sapply(selection_freq$all_sel, function(var) {
        parts <- strsplit(var, ":", fixed = TRUE)[[1]]
        cleaned <- sapply(parts, function(p) {
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
    }, USE.NAMES = FALSE)
    
    selection_freq <- distinct(selection_freq, all_sel, .keep_all = TRUE)
    
    cat("Completion rate:", n_boot_real / n_boot, "\n")
    
    list(
        selection_freq = selection_freq,
        selected_counts = selected_counts,
        n_boot_real = n_boot_real
    )
}

# Function to create and save variable frequency plot
plot_variable_frequency <- function(boot_forward_list, response_var, data_df, n_boot, output_dir = file.path(project_dir, "output/bootstrapped_models")) {
    results <- process_bootstrap_results(boot_forward_list, data_df, n_boot)
    
    p1 <- ggplot(
        results$selection_freq %>%
            rowwise() %>%
            mutate(all_sel = get_pretty_name(all_sel)) %>%
            ungroup() %>%
            head(30) %>%
            mutate(all_sel = factor(all_sel, levels = rev(all_sel))),
        aes(x = all_sel, y = (Freq / results$n_boot_real))
    ) +
        geom_bar(stat = "identity", fill = "steelblue") +
        xlab("Chosen Variable") +
        ylab(NULL) +
        hrbrthemes::scale_y_percent(limits = c(0, 1)) +
        ggtitle("Bootstrapped models - Variable Frequency") +
        theme_bw(12) +
        coord_flip()
    
    ggsave(p1, filename = paste0(output_dir, "/bootstrapped_models_variable_inclusion_", response_var, ".png"), width = 10, height = 6)
    
    return(results)
}

# Function to analyze and plot model combinations
plot_model_combinations <- function(boot_forward_list, response_var, data_df, n_boot, output_dir = file.path(project_dir, "output/bootstrapped_models")) {
    results <- process_bootstrap_results(boot_forward_list, data_df, n_boot)
    
    # Clean selected variables
    clean_selected_vars <- function(vars, model_vars) {
        sapply(vars, function(var) {
            if (grepl(":", var, fixed = TRUE)) {
                parts <- strsplit(var, ":", fixed = TRUE)[[1]]
                cleaned <- sapply(parts, function(p) {
                    hits <- agrep(p, names(model_vars), value = TRUE, max.distance = 20)
                    if (length(hits) > 0) {
                        d <- adist(p, hits)
                        best <- hits[which.min(d)][1]
                    } else {
                        best <- p
                    }
                    clean_var_names(best)
                })
                paste(sort(cleaned), collapse = ":")
            } else {
                hits <- agrep(var, names(model_vars), value = TRUE, max.distance = 20)
                if (length(hits) > 0) {
                    d <- adist(var, hits)
                    best <- hits[which.min(d)][1]
                } else {
                    best <- var
                }
                clean_var_names(best)
            }
        }, USE.NAMES = FALSE) %>% unique()
    }
    
    selected_counts_cleaned <- mclapply(results$selected_counts, function(vars) {
        if (is.null(vars)) {
            return(NULL)
        }
        clean_selected_vars(vars, data_df) %>% sort()
    }, mc.cores = detectCores() - 1)
    
    # Create binary matrix
    all_variables <- unique(unlist(selected_counts_cleaned)) %>% sort()
    binary_matrix <- sapply(selected_counts_cleaned, function(vars) {
        as.integer(all_variables %in% vars)
    })
    
    # Count combinations
    combination_counts <- table(apply(binary_matrix, 2, paste, collapse = ""))
    combination_counts <- sort(combination_counts, decreasing = TRUE)
    
    print(combination_counts)
    
    # Plot top combinations
    top_n <- max(5, sum(combination_counts > n_boot * 0.1))
    top_combos <- head(names(combination_counts), top_n)
    top_combo_counts <- combination_counts[top_combos]
    
    combo_matrix <- sapply(top_combos, function(combo) {
        as.logical(as.numeric(strsplit(combo, "")[[1]]))
    })
    rownames(combo_matrix) <- all_variables
    colnames(combo_matrix) <- paste0("Model ", seq_along(top_combos), " (", as.integer(top_combo_counts), "x)")
    
    
    
    combo_df <- as.data.frame(combo_matrix)
    combo_df$Variable <- rownames(combo_matrix)
    combo_long <- tidyr::pivot_longer(combo_df, -Variable, names_to = "Model", values_to = "Present")
    
    p <- ggplot(combo_long, aes(x = Model, y = Variable, fill = Present)) +
        geom_tile(color = "white") +
        scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey90")) +
        labs(
            title = paste("Top", top_n, "Model Combinations (by frequency)"),
            x = "Model Combination (count)",
            y = "Variable"
        ) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10)
        )
    
    print(p)
    
    ggsave(p, filename = paste0(output_dir, "/bootstrapped_models_combinations_", response_var, ".png"), width = 10, height = 6)
    
    return(list(
        combination_counts = combination_counts,
        all_variables = all_variables
    ))
}