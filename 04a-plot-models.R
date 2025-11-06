output_dir <- file.path(project_dir, "output/model_plots/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

get_pretty_name_v2 <- function(var) {
    # Remove .L, .C, .Q suffixes (for ordered factors)
    var_clean <- sub("\\.(L|C|Q)$", "", var)
    # Handle interaction terms like "Var1:Var2"
    if (grepl(":", var_clean, fixed = TRUE)) {
        parts <- strsplit(var_clean, ":", fixed = TRUE)[[1]]
        pretty_parts <- vapply(parts, get_pretty_name_v2, character(1))
        # Find and re-append any .L/.C/.Q suffixes from original var
        suffixes <- regmatches(var, gregexpr("\\.(L|C|Q)", var))[[1]]
        if (length(suffixes) == length(pretty_parts)) {
            pretty_parts <- paste0(pretty_parts, suffixes)
        }
        return(paste(pretty_parts, collapse = ":"))
    }
    # Exact match
    if (var_clean %in% pretty_names$variable) {
        pretty <- pretty_names$pretty_name[match(var_clean, pretty_names$variable)]
    } else {
        # Try to match a prefix in pretty_names$variable
        candidates <- pretty_names$variable[startsWith(var_clean, pretty_names$variable)]
        if (length(candidates) > 0) {
            # pick the longest matching prefix
            best <- candidates[which.max(nchar(candidates))]
            pretty_prefix <- pretty_names$pretty_name[pretty_names$variable == best]
            suffix <- substring(var_clean, nchar(best) + 1)
            # strip leading separators from suffix
            suffix <- sub("^[ _:\\.]+", "", suffix)
            pretty <- paste0(pretty_prefix, if (nzchar(suffix)) {
                paste0(" (", suffix, ")")
            } else {
                ""
            })
        } else {
            pretty <- var_clean
        }
    }
    # Re-append .L/.C/.Q if present in original var
    suffix <- sub(".*(\\.(L|C|Q))$", "\\1", var)
    if (grepl("\\.(L|C|Q)$", var)) {
        # Map suffix to label
        suffix_label <- switch(suffix,
            ".L" = " [linear]",
            ".C" = " [cubic]",
            ".Q" = " [quadratic]",
            ""
        )
        pretty <- paste0(pretty, suffix_label)
    }
    pretty
}

# Helper function to plot model coefficients with pretty names, AI/Path coloring, and R² box
plot_model_coefficients <- function(model, importance_summary, get_pretty_name_v2, title = "TCC Discrepancy Model") {
    coef_names <- names(model$coefficients[-1])
    orig_vars <- importance_summary$orig_vars

    # Remove .L/.Q/.C suffixes for matching
    coef_base <- sub("\\.(L|C|Q)$", "", coef_names)

    # Remove factor level suffixes like "High" or "Good" from coef_base and orig_vars for matching
    remove_level_suffix <- function(x) {
        # Remove suffix if it matches any in the list (e.g., "High", "Good") at the end of the string
        suffixes <- c("High", "Good")
        pattern <- paste0("(", paste(suffixes, collapse = "|"), ")$")
        sub(pattern, "", x)
    }
    coef_base <- remove_level_suffix(coef_base)



    # For each orig_var, find matching coef indices (with or without .L/.Q)
    ordered_indices <- unlist(lapply(orig_vars, function(var) which(coef_base == var)))
    ordered_indices <- ordered_indices[ordered_indices > 0]

    # Compute y-axis limits so that min and max are symmetric around zero
    coef_vals <- model$coefficients
    coef_vals <- coef_vals[!is.na(coef_vals)]
    coef_names <- coef_names[which(!is.na(coef_vals))-1] # -1 for intercept already removed
    # Get standard errors for coefficients
    coef_se <- summary(model)$coefficients[, "Std. Error"]
    # Compute range including SE
    coef_range <- range(coef_vals + coef_se, coef_vals - coef_se, na.rm = TRUE)
    abs_max <- max(abs(coef_range))
    y_limits <- c(-abs_max - 5, abs_max + 5)

    pm <- plot_model(
        model,
        show.intercept = TRUE, order.terms = c(1, (ordered_indices + 1)),
        axis.labels = c("(Intercept)" = "Intercept", vapply(
            coef_names,
            get_pretty_name_v2,
            character(1)
        ))
    ) + theme_bw(18) + ylab("Coefficient Estimate") + ggtitle(title) +
        scale_y_continuous(limits = y_limits)

    pm <- pm + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")

    # Annotator group: Intercept, AI, Pathologist
    pm$data$group <- c(
        "Intercept",
        grepl("^\\(?AI", vapply(coef_names, get_pretty_name_v2, character(1)))
    )

    scale_intercept_ai_path <- c("Intercept" = "#797979", "FALSE" = "#A767FF", "TRUE" = "#1041FF")

    pm <- pm +
        scale_color_manual(values = scale_intercept_ai_path) +
        scale_fill_manual(values = scale_intercept_ai_path)

    pm$layers[[2]]$aes_params$size <- 4

    pm <- pm +
        guides(color = guide_legend(title = "Annotator", override.aes = list(size = 4))) +
        scale_color_manual(
            values = scale_intercept_ai_path[c("TRUE", "FALSE")],
            breaks = c("TRUE", "FALSE"),
            labels = c("AI", "Pathologist")
        ) +
        theme(
            legend.position = "bottom"
        )

    # Calculate R² for the model
    r2 <- summary(model)$r.squared
    r2 <- round(r2 * 100, 0) # Convert to percentage

    # Add R² box to the pm plot
    # Extract y-axis limits from the plot
    gb <- ggplot_build(pm)
    y_range <- gb$layout$panel_params[[1]]$x.range # sic
    # Round ymax down to the nearest multiple of 5
    ymax <- y_range[2] - (y_range[2] / 5)
    ymin <- y_range[1]

    # Box size as fraction of y-axis range
    box_height <- (ymax - ymin) * 0.18
    box_width <- 1 # width in terms of x-axis (1 coefficient bar width)
    # Place box at top right (right of first coefficient bar)
    xmin <- 2
    xmax <- xmin + box_width
    box_ymin <- ymax - box_height
    box_ymax <- ymax

    pm <- pm +
        annotate(
            "rect",
            xmin = xmin, xmax = xmax, ymin = box_ymin, ymax = box_ymax,
            fill = "#72c2ff", alpha = 0.85, color = "black"
        ) +
        annotate(
            "text",
            x = (xmin + xmax) / 2, y = (box_ymin + box_ymax) / 2,
            label = paste0(sprintf("R² = %.0f", r2), "%"),
            size = 5, fontface = "bold", color = "black"
        )

    pm
}

# Example usage for Path vs AI model
pm <- plot_model_coefficients(
    classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$model,
    classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$importance_summary,
    get_pretty_name_v2,
    title = "TCC Discrepancy Pathologist vs AI"
)
ggsave(paste0(output_dir, "TCC_discrepancy_model.pdf"), pm, width = 8, height = 10, dpi = 300)

# Example usage for FMI vs AI model
pm_fmi_vs_ai <- plot_model_coefficients(
    classic_no_forced_interactions$TCC_FMI_minus_TCC_AI$model,
    classic_no_forced_interactions$TCC_FMI_minus_TCC_AI$importance_summary,
    get_pretty_name_v2,
    title = "TCC Discrepancy FMI vs AI"
)
ggsave(paste0(output_dir, "TCC_discrepancy_model_fmi_vs_ai.pdf"), pm_fmi_vs_ai, width = 8, height = 9, dpi = 300)

# Example usage for Path vs FMI model
pm_path_vs_fmi <- plot_model_coefficients(
    classic_no_forced_interactions$TCC_Patho_minus_TCC_FMI$model,
    classic_no_forced_interactions$TCC_Patho_minus_TCC_FMI$importance_summary,
    get_pretty_name_v2,
    title = "TCC Discrepancy Pathologist vs FMI"
)

ggsave(paste0(output_dir, "TCC_discrepancy_model_path_vs_fmi.pdf"), pm_path_vs_fmi, width = 8, height = 8, dpi = 300)





# Helper function to get sensible grouping values for a variable
get_grouping_values <- function(var_name, data, n_groups = 4) {
    var_data <- data[[var_name]]
    if (is.factor(var_data)) {
        # For factors, return all levels (or a subset if too many)
        levs <- levels(var_data)
        if (length(levs) <= n_groups) return(levs)
        # Sample evenly across levels
        idx <- round(seq(1, length(levs), length.out = n_groups))
        return(levs[idx])
    } else {
        # For numeric, use quantiles excluding extremes
        quantiles <- quantile(var_data, probs = seq(0.1, 0.9, length.out = n_groups), na.rm = TRUE)
        # Round to sensible values
        range_val <- max(quantiles) - min(quantiles)
        if (range_val > 1000) {
            quantiles <- round(quantiles, -2) # Round to nearest 100
        } else if (range_val > 100) {
            quantiles <- round(quantiles, -1) # Round to nearest 10
        } else {
            quantiles <- round(quantiles, 1)
        }
        return(unique(quantiles))
    }
}

# Loop through all models in classic_no_forced_interactions
for (model_name in names(classic_no_forced_interactions)) {
    cat("Processing model:", model_name, "\n")
    
    model_obj <- classic_no_forced_interactions[[model_name]]$model
    model_data <- data_df_renamed
    
    # Get predictor names (excluding intercept)
    predictor_names <- names(model_obj$coefficients)[-1]
    
    # Remove .L, .Q, .C suffixes and factor level suffixes to get base variable names
    get_base_var <- function(pred_name) {
        # Remove polynomial contrast suffixes
        pred_clean <- sub("\\.(L|Q|C)$", "", pred_name)
        # Remove common factor level suffixes
        pred_clean <- sub("(High|Good|Low|Bad|Medium)$", "", pred_clean)
        return(pred_clean)
    }
    
    # Process each predictor
    processed_vars <- character(0)
    
    for (pred in predictor_names) {
        # Skip if already processed (handles .L, .Q, .C variants)
        base_pred <- get_base_var(pred)
        if (base_pred %in% processed_vars) next
        
        # Check if interaction term
        if (grepl(":", pred, fixed = TRUE)) {
            # Interaction term
            parts <- strsplit(base_pred, ":", fixed = TRUE)[[1]]
            if (length(parts) != 2) next
            
            var_base <- parts[1]
            var_base_2 <- parts[2]
            
            # Get appropriate grouping values for var_base_2
            var_base_2_values <- get_grouping_values(var_base_2, model_data, n_groups = 4)
            
            tryCatch({
                p <- plot_main_effect(
                    model = model_obj,
                    model_data = model_data,
                    var_base = var_base,
                    var_base_2 = var_base_2,
                    var_base_2_values = var_base_2_values,
                    get_pretty_name_v2 = get_pretty_name_v2,
                    unscale_center = TRUE
                ) + 
                    theme(legend.position = "bottom") + 
                    ggtitle(model_name) +
                    ylab("Predicted Discrepancy")
                
                # Create safe filename
                safe_filename <- paste0(
                    output_dir,
                    model_name, "_interaction_",
                    make.names(var_base), "_by_",
                    make.names(var_base_2), ".pdf"
                )
                
                ggsave(safe_filename, p, width = 10, height = 8, dpi = 300)
                cat("  Saved interaction plot:", safe_filename, "\n")
            }, error = function(e) {
                cat("  Error plotting interaction", var_base, ":", var_base_2, "-", e$message, "\n")
            })
            
            processed_vars <- c(processed_vars, base_pred)
            
        } else {
            # Main effect term
            var_base <- base_pred
            
            tryCatch({
                p <- plot_main_effect(
                    model = model_obj,
                    model_data = model_data,
                    var_base = var_base,
                    get_pretty_name_v2 = get_pretty_name_v2,
                    unscale_center = TRUE
                ) + 
                    ggtitle(model_name) +
                    ylab("Predicted Discrepancy")
                
                # Create safe filename
                safe_filename <- paste0(
                    output_dir,
                    model_name, "_main_effect_",
                    make.names(var_base), ".pdf"
                )
                
                ggsave(safe_filename, p, width = 8, height = 8, dpi = 300)
                cat("  Saved main effect plot:", safe_filename, "\n")
            }, error = function(e) {
                cat("  Error plotting main effect", var_base, "-", e$message, "\n")
            })
            
            processed_vars <- c(processed_vars, var_base)
        }
    }
}

cat("All plots generated successfully!\n")

