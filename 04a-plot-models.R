output_dir <- "discrepancies/2025-10-Data-Version/model_plots/"
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
            pretty <- paste0(pretty_prefix, if (nzchar(suffix)) suffix else "")
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
ggsave(paste0(output_dir, "TCC_discrepancy_model.pdf"), pm, width = 8, height = 8, dpi = 300)

# Example usage for FMI vs AI model
pm_fmi_vs_ai <- plot_model_coefficients(
    classic_no_forced_interactions$TCC_FMI_minus_TCC_AI$model,
    classic_no_forced_interactions$TCC_FMI_minus_TCC_AI$importance_summary,
    get_pretty_name_v2,
    title = "TCC Discrepancy FMI vs AI"
)
ggsave(paste0(output_dir, "TCC_discrepancy_model_fmi_vs_ai.pdf"), pm_fmi_vs_ai, width = 8, height = 8, dpi = 300)

# Example usage for Path vs FMI model
pm_path_vs_fmi <- plot_model_coefficients(
    classic_no_forced_interactions$TCC_Patho_minus_TCC_FMI$model,
    classic_no_forced_interactions$TCC_Patho_minus_TCC_FMI$importance_summary,
    get_pretty_name_v2,
    title = "TCC Discrepancy Pathologist vs FMI"
)

ggsave(paste0(output_dir, "TCC_discrepancy_model_path_vs_fmi.pdf"), pm_path_vs_fmi, width = 8, height = 8, dpi = 300)

# main effects
# Function to plot main effect for a variable, with optional grouping by a second variable
plot_main_effect <- function(
    model, model_data, var_base, get_pretty_name_v2 = identity, ylab = "Predicted Discrepancy",
    title_prefix = "Main Effect of", var_base_2 = NULL, var_base_2_values = NULL, palette = NULL,
    var_base_2_correction = NULL, var_base_2_correction_values = NULL,
    unscale_center = FALSE, y_limits = c(-100, 100)) {
    if (unscale_center) {
        raw_data <- data_df_pre_scaling
        names(raw_data) <- make.names(names(raw_data), unique = TRUE)
    }

    # Small helpers to map raw -> model scale (y_model ~ a + b * x_raw)
    fit_lin_scaler <- function(model_col, raw_col) {
        idx <- is.finite(model_col) & is.finite(raw_col)
        if (!any(idx) || sum(idx) < 2) {
            return(list(a = 0, b = 1))
        }
        co <- coef(lm(model_col[idx] ~ raw_col[idx]))
        list(a = unname(co[[1]]), b = unname(co[[2]]))
    }
    apply_scaler <- function(raw_vals, scaler) scaler$a + scaler$b * raw_vals
    raw_log <- function(x) ifelse(x > 1, log(x) + 1, x)

    model_vars <- names(model$model)

    # Color-blind friendly palette (Okabe-Ito)
    cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

    # NEW: support factor var_base
    var_is_factor <- is.factor(model_data[[var_base]])
    print(head(model_data[[var_base]]))
    print(var_is_factor)
    if (var_is_factor) {
        if (!(var_base %in% model_vars)) stop("No matching variable found in model for ", var_base)

        levels_base <- levels(model$model[[var_base]])
        # Prepare combinations across var_base levels and (optional) var_base_2 values
        if (is.null(var_base_2)) {
            combos <- tibble::tibble(level = levels_base)
            var2_is_factor <- FALSE
        } else {
            if (is.null(var_base_2_values)) {
                stop("If var_base_2 is provided, var_base_2_values must also be provided.")
            }
            var2_is_factor <- is.factor(model$model[[var_base_2]])
            if (var2_is_factor) {
                # ensure provided values are valid factor levels
                valid_levels <- levels(model$model[[var_base_2]])
                if (!all(var_base_2_values %in% valid_levels)) {
                    stop("var_base_2_values must be a subset of factor levels of var_base_2")
                }
            }
            combos <- tidyr::expand_grid(
                level = levels_base,
                val2 = var_base_2_values
            )
        }

        # Pre-fit scalers for numeric var_base_2 (only if unscale_center)
        if (!is.null(var_base_2) && !var2_is_factor && unscale_center) {
            s_base2 <- fit_lin_scaler(model_data[[var_base_2]], raw_data[[var_base_2]])
            if (paste0("log_", var_base_2) %in% model_vars) {
                s_log_base2 <- fit_lin_scaler(model_data[[paste0("log_", var_base_2)]], raw_log(raw_data[[var_base_2]]))
            }
            if (paste0("sq_", var_base_2) %in% model_vars) {
                s_sq_base2 <- fit_lin_scaler(model_data[[paste0("sq_", var_base_2)]], raw_data[[var_base_2]]^2)
            }
        }

        # Helper: build a 1-row newdata at defaults, then set var_base and var_base_2 as required
        build_newdata <- function(level_val, var2_val = NULL) {
            nd <- list()
            for (v in model_vars) {
                if (v == var_base) {
                    nd[[v]] <- factor(level_val, levels = levels(model$model[[v]]))
                } else if (!is.null(var_base_2) && v == var_base_2) {
                    if (var2_is_factor) {
                        nd[[v]] <- factor(var2_val, levels = levels(model$model[[v]]))
                    } else {
                        if (unscale_center) {
                            nd[[v]] <- apply_scaler(var2_val, s_base2)
                        } else {
                            nd[[v]] <- var2_val
                        }
                    }
                } else if (is.factor(model$model[[v]])) {
                    nd[[v]] <- factor(levels(model$model[[v]])[1], levels = levels(model$model[[v]]))
                } else {
                    nd[[v]] <- stats::median(model_data[[v]], na.rm = TRUE)
                }
            }
            # Fill log_/sq_ for numeric var_base_2 if present
            if (!is.null(var_base_2) && !var2_is_factor) {
                if (paste0("log_", var_base_2) %in% model_vars) {
                    nd[[paste0("log_", var_base_2)]] <- if (unscale_center) {
                        apply_scaler(raw_log(var2_val), s_log_base2)
                    } else {
                        ifelse(var2_val > 1, log(var2_val) + 1, var2_val)
                    }
                }
                if (paste0("sq_", var_base_2) %in% model_vars) {
                    nd[[paste0("sq_", var_base_2)]] <- if (unscale_center) {
                        apply_scaler(var2_val^2, s_sq_base2)
                    } else {
                        var2_val^2
                    }
                }
            }
            as.data.frame(nd, stringsAsFactors = FALSE)
        }

        # Simulate draws from the predictive uncertainty (same basis as numeric ribbon: se.fit)
        draw_n <- 1000
        preds <- purrr::pmap_dfr(
            combos,
            function(level, val2 = NULL) {
                nd <- build_newdata(level, val2)
                pr <- predict(model, nd, se.fit = TRUE)
                tibble::tibble(
                    var_level = level,
                    group = if (is.null(var_base_2)) NA_character_ else as.character(val2),
                    pred = stats::rnorm(draw_n, mean = pr$fit[1], sd = pr$se.fit[1])
                )
            }
        )

        preds$var_level <- factor(preds$var_level, levels = levels_base)

        # Build plot: boxplots per level (dodged by var_base_2 when provided)
        p <- ggplot2::ggplot(
            preds,
            ggplot2::aes(
                x = var_level,
                y = pred,
                fill = if (!is.null(var_base_2)) group
            )
        ) +
            ggplot2::geom_boxplot(
                outlier.shape = NA, alpha = 0.7,
                position = if (is.null(var_base_2)) "identity" else ggplot2::position_dodge2(width = 0.8, preserve = "single")
            ) +
            ggplot2::labs(
                x = get_pretty_name_v2(var_base),
                y = ylab,
                title = paste(
                    title_prefix,
                    get_pretty_name_v2(var_base),
                    if (!is.null(var_base_2)) paste("by", get_pretty_name_v2(var_base_2)) else ""
                )
            ) +
            ggplot2::theme_bw(14) +
            ggplot2::coord_cartesian(ylim = y_limits)

        if (is.null(var_base_2)) {
            p <- p + ggplot2::guides(fill = "none")
        } else {
            pal_vals <- if (is.null(palette)) cb_palette[seq_along(var_base_2_values)] else palette
            p <- p +
                ggplot2::scale_fill_manual(
                    values = pal_vals,
                    name = get_pretty_name_v2(var_base_2),
                    labels = if (var2_is_factor) var_base_2_values else scales::comma(var_base_2_values)
                )
        }

        return(p)
    }

    # --- Numeric var_base path (existing code) ---
    # Find all columns in the model that are the base variable, log_ or sq_ transformed
    var_names <- c(var_base, paste0("log_", var_base), paste0("sq_", var_base))
    present_vars <- var_names[var_names %in% model_vars]
    if (length(present_vars) == 0) stop("No matching variable found in model for ", var_base)

    # Get the range of the base variable (raw if requested)
    var_range <- if (unscale_center) range(raw_data[[var_base]], na.rm = TRUE) else range(model_data[[var_base]], na.rm = TRUE)
    n_points <- 1000

    # If no grouping variable
    if (is.null(var_base_2)) {
        # Build newdata with correct model-scale values; keep x_raw for plotting
        x_raw <- if (unscale_center) seq(var_range[1], var_range[2], length.out = n_points) else NULL
        newdata <- data.frame(matrix(ncol = 1, nrow = n_points))
        colnames(newdata) <- var_base

        # Pre-fit scalers if needed
        if (unscale_center) {
            s_base <- fit_lin_scaler(model_data[[var_base]], raw_data[[var_base]])
            # log term scaler over dataset
            if (paste0("log_", var_base) %in% model_vars) {
                s_log_base <- fit_lin_scaler(model_data[[paste0("log_", var_base)]], raw_log(raw_data[[var_base]]))
            }
            if (paste0("sq_", var_base) %in% model_vars) {
                s_sq_base <- fit_lin_scaler(model_data[[paste0("sq_", var_base)]], raw_data[[var_base]]^2)
            }
        }

        for (v in model_vars) {
            if (v == var_base) {
                if (unscale_center) {
                    newdata[[v]] <- apply_scaler(x_raw, s_base)
                } else {
                    newdata[[v]] <- seq(var_range[1], var_range[2], length.out = n_points)
                }
            } else if (is.factor(model$model[[v]])) {
                newdata[[v]] <- factor(rep(levels(model$model[[v]])[1], n_points),
                    levels = levels(model$model[[v]])
                )
            } else {
                newdata[[v]] <- median(model_data[[v]], na.rm = TRUE)
            }
        }

        # Fill in log_ and sq_ columns if present based on raw x, then scale to model
        if (paste0("log_", var_base) %in% model_vars) {
            if (unscale_center) {
                newdata[[paste0("log_", var_base)]] <- apply_scaler(raw_log(x_raw), s_log_base)
            } else {
                newdata[[paste0("log_", var_base)]] <- ifelse(newdata[[var_base]] > 1, log(newdata[[var_base]]) + 1, newdata[[var_base]])
            }
        }
        if (paste0("sq_", var_base) %in% model_vars) {
            if (unscale_center) {
                newdata[[paste0("sq_", var_base)]] <- apply_scaler(x_raw^2, s_sq_base)
            } else {
                newdata[[paste0("sq_", var_base)]] <- newdata[[var_base]]^2
            }
        }

        pred <- predict(model, newdata, se.fit = TRUE)
        pred_df <- tibble(
            x = if (unscale_center) x_raw else newdata[[var_base]],
            fit = pred$fit,
            lwr = pred$fit - 1.96 * pred$se.fit,
            upr = pred$fit + 1.96 * pred$se.fit
        )

        p <- ggplot(pred_df, aes(x = x, y = fit)) +
            geom_line(color = "black", size = 1.2) +
            geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "black", alpha = 0.2) +
            labs(
                x = get_pretty_name_v2(var_base),
                y = ylab,
                title = paste(title_prefix, get_pretty_name_v2(var_base))
            ) +
            scale_x_continuous(labels = scales::comma) +
            theme_bw(14) +
            coord_cartesian(ylim = y_limits)
        return(p)
    }

    # If grouping variable is provided
    if (is.null(var_base_2_values)) {
        stop("If var_base_2 is provided, var_base_2_values must also be provided.")
    }

    var2_names <- c(var_base_2, paste0("log_", var_base_2), paste0("sq_", var_base_2))
    present_vars2 <- var2_names[var2_names %in% model_vars]
    if (length(present_vars2) == 0) stop("No matching variable found in model for ", var_base_2)

    if (is.null(palette)) palette <- cb_palette[seq_along(var_base_2_values)]

    if (!is.null(var_base_2_correction)) {
        if (is.null(var_base_2_correction_values)) {
            stop("If var_base_2_correction is provided, var_base_2_correction_values must also be provided.")
        }
        if (length(var_base_2_correction_values) != length(var_base_2_values)) {
            stop("var_base_2_correction_values must be the same length as var_base_2_values.")
        }
        var2_corr_names <- c(var_base_2_correction, paste0("log_", var_base_2_correction), paste0("sq_", var_base_2_correction))
        present_vars2_corr <- var2_corr_names[var2_names %in% model_vars]
        if (length(present_vars2_corr) == 0) stop("No matching variable found in model for ", var_base_2_correction)
    }

    # Pre-fit scalers for base variables if needed
    if (unscale_center) {
        s_base <- fit_lin_scaler(model_data[[var_base]], raw_data[[var_base]])
        if (paste0("log_", var_base) %in% model_vars) {
            s_log_base <- fit_lin_scaler(model_data[[paste0("log_", var_base)]], raw_log(raw_data[[var_base]]))
        }
        if (paste0("sq_", var_base) %in% model_vars) {
            s_sq_base <- fit_lin_scaler(model_data[[paste0("sq_", var_base)]], raw_data[[var_base]]^2)
        }
        s_base2 <- fit_lin_scaler(model_data[[var_base_2]], raw_data[[var_base_2]])
        if (paste0("log_", var_base_2) %in% model_vars) {
            s_log_base2 <- fit_lin_scaler(model_data[[paste0("log_", var_base_2)]], raw_log(raw_data[[var_base_2]]))
        }
        if (paste0("sq_", var_base_2) %in% model_vars) {
            s_sq_base2 <- fit_lin_scaler(model_data[[paste0("sq_", var_base_2)]], raw_data[[var_base_2]]^2)
        }
        if (!is.null(var_base_2_correction)) {
            s_corr <- fit_lin_scaler(model_data[[var_base_2_correction]], raw_data[[var_base_2_correction]])
            if (paste0("log_", var_base_2_correction) %in% model_vars) {
                s_log_corr <- fit_lin_scaler(model_data[[paste0("log_", var_base_2_correction)]], raw_log(raw_data[[var_base_2_correction]]))
            }
            if (paste0("sq_", var_base_2_correction) %in% model_vars) {
                s_sq_corr <- fit_lin_scaler(model_data[[paste0("sq_", var_base_2_correction)]], raw_data[[var_base_2_correction]]^2)
            }
        }
    }

    all_pred <- lapply(seq_along(var_base_2_values), function(i) {
        val2_raw <- var_base_2_values[i]
        x_raw <- if (unscale_center) seq(var_range[1], var_range[2], length.out = n_points) else NULL

        newdata <- data.frame(matrix(ncol = 1, nrow = n_points))
        colnames(newdata) <- var_base

        for (v in model_vars) {
            if (v == var_base) {
                if (unscale_center) {
                    newdata[[v]] <- apply_scaler(x_raw, s_base)
                } else {
                    newdata[[v]] <- seq(var_range[1], var_range[2], length.out = n_points)
                }
            } else if (v == var_base_2) {
                if (unscale_center) {
                    newdata[[v]] <- rep(apply_scaler(val2_raw, s_base2), n_points)
                } else {
                    newdata[[v]] <- rep(val2_raw, n_points)
                }
            } else if (!is.null(var_base_2_correction) && v == var_base_2_correction) {
                if (unscale_center) {
                    newdata[[v]] <- rep(apply_scaler(var_base_2_correction_values[i], s_corr), n_points)
                } else {
                    newdata[[v]] <- rep(var_base_2_correction_values[i], n_points)
                }
            } else if (is.factor(model$model[[v]])) {
                newdata[[v]] <- factor(rep(levels(model$model[[v]])[1], n_points),
                    levels = levels(model$model[[v]])
                )
            } else {
                newdata[[v]] <- median(model_data[[v]], na.rm = TRUE)
            }
        }

        # log/sq for var_base from raw x
        if (paste0("log_", var_base) %in% model_vars) {
            if (unscale_center) {
                newdata[[paste0("log_", var_base)]] <- apply_scaler(raw_log(x_raw), s_log_base)
            } else {
                newdata[[paste0("log_", var_base)]] <- ifelse(newdata[[var_base]] > 1, log(newdata[[var_base]]) + 1, newdata[[var_base]])
            }
        }
        if (paste0("sq_", var_base) %in% model_vars) {
            if (unscale_center) {
                newdata[[paste0("sq_", var_base)]] <- apply_scaler(x_raw^2, s_sq_base)
            } else {
                newdata[[paste0("sq_", var_base)]] <- newdata[[var_base]]^2
            }
        }

        # log/sq for var_base_2 from raw val2
        if (paste0("log_", var_base_2) %in% model_vars) {
            if (unscale_center) {
                newdata[[paste0("log_", var_base_2)]] <- rep(apply_scaler(raw_log(val2_raw), s_log_base2), n_points)
            } else {
                newdata[[paste0("log_", var_base_2)]] <- rep(ifelse(val2_raw > 1, log(val2_raw) + 1, val2_raw), n_points)
            }
        }
        if (paste0("sq_", var_base_2) %in% model_vars) {
            if (unscale_center) {
                newdata[[paste0("sq_", var_base_2)]] <- rep(apply_scaler(val2_raw^2, s_sq_base2), n_points)
            } else {
                newdata[[paste0("sq_", var_base_2)]] <- rep(val2_raw^2, n_points)
            }
        }

        # log/sq for correction var
        if (!is.null(var_base_2_correction)) {
            val2_corr_raw <- var_base_2_correction_values[i]
            if (paste0("log_", var_base_2_correction) %in% model_vars) {
                if (unscale_center) {
                    newdata[[paste0("log_", var_base_2_correction)]] <- rep(apply_scaler(raw_log(val2_corr_raw), s_log_corr), n_points)
                } else {
                    newdata[[paste0("log_", var_base_2_correction)]] <- rep(ifelse(val2_corr_raw > 1, log(val2_corr_raw) + 1, val2_corr_raw), n_points)
                }
            }
            if (paste0("sq_", var_base_2_correction) %in% model_vars) {
                if (unscale_center) {
                    newdata[[paste0("sq_", var_base_2_correction)]] <- rep(apply_scaler(val2_corr_raw^2, s_sq_corr), n_points)
                } else {
                    newdata[[paste0("sq_", var_base_2_correction)]] <- rep(val2_corr_raw^2, n_points)
                }
            }
        }

        pred <- predict(model, newdata, se.fit = TRUE)
        tibble(
            x = if (unscale_center) x_raw else newdata[[var_base]],
            fit = pred$fit,
            lwr = pred$fit - 1.96 * pred$se.fit,
            upr = pred$fit + 1.96 * pred$se.fit,
            group = as.character(val2_raw)
        )
    })
    pred_df <- bind_rows(all_pred)

    p <- ggplot(pred_df, aes(x = x, y = fit, color = group, fill = group)) +
        geom_line(size = 1.2) +
        geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
        scale_color_manual(
            values = palette,
            name = get_pretty_name_v2(var_base_2),
            labels = scales::comma(var_base_2_values)
        ) +
        scale_fill_manual(
            values = palette,
            name = get_pretty_name_v2(var_base_2),
            labels = scales::comma(var_base_2_values)
        ) +
        scale_x_continuous(labels = scales::comma) +
        labs(
            x = get_pretty_name_v2(var_base),
            y = ylab,
            title = paste(title_prefix, get_pretty_name_v2(var_base), "by", get_pretty_name_v2(var_base_2))
        ) +
        theme_bw(14) +
        coord_cartesian(ylim = c(-100, 100))
    return(p)
}
# Example usage:
# Single variable:
me2 <- plot_main_effect(
    model = classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$model,
    model_data = data_df_renamed,
    var_base = "AI_total_necrosis_area_.",
    get_pretty_name_v2, unscale_center = TRUE
) + ggtitle("... vs AI") +
    ylab("Predicted Discrepancy")

# Function to add a value box in the bottom right corner of a ggplot
add_importance_box <- function(plot, imp_val, label = "Importance") {
    # Get importance value (in %)
    imp_val <- round(imp_val * 100, 0)
    # Get axis limits from plot
    gb <- ggplot_build(plot)
    xlim <- gb$layout$panel_params[[1]]$x.range
    ylim <- gb$layout$panel_params[[1]]$y.range
    # Box size as fraction of axis range
    box_width <- diff(xlim) * 0.13
    box_height <- diff(ylim) * 0.13
    # Box position (bottom right)
    xmin <- xlim[2] - box_width
    xmax <- xlim[2]
    ymin <- ylim[1]
    ymax <- ylim[1] + box_height
    # Add box and text
    plot +
        annotate("rect",
            xmin = xmin, xmax = xmax, ymin = ymin + 5, ymax = ymax,
            fill = "#f7e6ff", color = "black", alpha = 0.85
        ) +
        annotate("text",
            x = (xmin + xmax) / 2, y = (ymin + 5 + ymax) / 2,
            label = paste0(label, imp_val, "%"),
            size = 5, fontface = "bold", color = "black", vjust = 0.5
        )
}

# Example usage for me2:
# me2 <- add_importance_box(
#     me2,
#     imp_val=classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$importance_summary$RelativeImportance[classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$importance_summary$orig_vars=="log_TumorDetect.v1.2.0...Supporting.Result...Area.of.Necrosis...mm.2"],
#     label = ""
# )

# With grouping variable:
me1 <- plot_main_effect(
    model = classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$model,
    model_data = data_df_renamed,
    var_base = "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cells.in.normal.tissue",
    var_base_2 = "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cancer.cells.in.tumor",
    var_base_2_values = c(1000, 10000, 100000, 1000000),
    get_pretty_name_v2,
    unscale_center = TRUE
) + theme(legend.position = "bottom") + ggtitle("... vs AI") +
    ylab("Predicted Discrepancy")


# Add importance box to me1
# me1 <- add_importance_box(
#     me1,
#     imp_val=classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$importance_summary$RelativeImportance[classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$importance_summary$orig_vars=="Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cells.in.normal.tissue"],
#     label = ""
# )

# vs fmi model effects


# Single variable:
me4 <- plot_main_effect(
    model = path_vs_fmi$model,
    model_data = all_model_variables_renamed %>% mutate(TumorDetect.v1.2.0...Supporting.Result...Area.of.Necrosis...mm.2 = data$`TumorDetect v1.2.0 - Supporting Result - Area of Necrosis - mm^2`),
    var_base = "TumorDetect.v1.2.0...Supporting.Result...Area.of.Necrosis...mm.2",
    get_pretty_name_v2
) + ggtitle("... vs FMI") +
    ylab("Predicted Discrepancy")

# With grouping variable:
me3 <- ggplot() +
    theme_void()
# outdated and not updated because not needed
# plot_main_effect(
#     model = path_vs_fmi$model,
#     model_data = all_model_variables_renamed,
#     var_base = "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cells.in.normal.tissue",
#     var_base_2 = "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cancer.cells.in.tumor",
#     var_base_2_values = c(1000, 10000, 100000, 1000000),
#     var_base_2_correction = "log_Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.area.of.cancer...mm.2",
#     var_base_2_correction_values = log(c(1000+1, 10000+1, 100000+1, 1000000+1))*0.5,
#     get_pretty_name
# ) + theme(legend.position = "bottom") + ggtitle("... vs FMI") +
#     ylab("Predicted Discrepancy")


# Arrange plots: pm on the left 50%, me1/me2 top right, me3/me4 bottom right
library(patchwork)

# Remove legends from individual plots
me1_noleg <- me1 + theme(legend.position = "none") + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")
me2_noleg <- me2 + theme(legend.position = "none") + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")
me3_noleg <- me3 + theme(legend.position = "none") + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")
me4_noleg <- me4 + theme(legend.position = "none") + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")

# Extract a legend from one of the grouped plots (me1 or me3)
get_legend <- function(p) {
    # Extract only the legend grob from a ggplot object
    g <- ggplotGrob(p + theme(legend.position = "bottom"))
    legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
    if (length(legend_index) == 0) stop("No legend found in plot")
    patchwork::wrap_elements(g$grobs[[legend_index]])
}
legend_me1 <- get_legend(me1)

legend_pm <- get_legend(pm)

pm_noleg <- pm + theme(legend.position = "none")

# Compose the main plot without legends

left_hand <- (pm_noleg / legend_pm + plot_layout(heights = c(0.9, 0.1)))
right_hand <- ((me1_noleg + me2_noleg) / (me3_noleg + me4_noleg) / legend_me1 + plot_layout(heights = c(0.45, 0.45, 0.1)))
main_plot <- (left_hand | right_hand) + plot_layout(widths = c(0.4, 0.6))
main_plot






#### other (?) code

## plot each main effect in a grid

me1 <- plot_main_effect(
    model = classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$model,
    model_data = all_model_variables_renamed,
    var_base = "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cells.in.normal.tissue",
    var_base_2 = "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cancer.cells.in.tumor",
    var_base_2_values = c(1000, 10000, 100000, 1000000),
    get_pretty_name_v2,
    unscale_center = TRUE
) + theme(legend.position = "bottom") + ggtitle("... vs AI") +
    ylab("Predicted Discrepancy")

# Get coefficient names, exclude intercept
coef_names <- names(classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$model$coefficients)[-1]

# Remove "log_" or "sq_" prefixes and get unique base variables
clean_names <- str_replace(coef_names, "^(log_|sq_)", "")
clean_names <- str_replace(clean_names, "\\.(L|C|Q)$", "")
unique_clean_names <- unique(clean_names)

# Exclude specific variables
exclude_vars <- c(
    "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cells.in.normal.tissue",
    "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cancer.cells.in.tumor"
)
final_vars <- setdiff(unique_clean_names, exclude_vars)

# Generate plots for each final variable
main_effect_plots <- lapply(final_vars, function(var) {
    plot_main_effect(
        model = classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$model,
        model_data = all_model_variables_renamed,
        var_base = var,
        get_pretty_name_v2,
        unscale_center = TRUE,
        y_limits = c(-25, 25)
    ) + ggtitle(NULL) +
        ylab(NULL)
})

# Optionally, name the list elements for reference
names(main_effect_plots) <- final_vars

# Save each main effect plot as a PNG file
for (i in seq_along(main_effect_plots)) {
    var_name <- names(main_effect_plots)[i]
    file_name <- paste0("test_", var_name, ".png")
    ggsave(filename = file_name, plot = main_effect_plots[[i]], width = 8, height = 6, dpi = 300)
}

me1 <- plot_main_effect(
    model = classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$model,
    model_data = all_model_variables_renamed,
    var_base = "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cells.in.normal.tissue",
    var_base_2 = "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cancer.cells.in.tumor",
    var_base_2_values = c(1000, 10000, 100000, 1000000),
    get_pretty_name_v2,
    unscale_center = TRUE
) + theme(legend.position = "bottom") + ggtitle(NULL) +
    ylab(NULL) +
    scale_y_continuous(breaks = seq(-100, 100, 20))

ggsave("test.png", me1, width = 8, height = 6, dpi = 300)

# Combine me1 on top and the four main_effect_plots in a 2x2 grid underneath
top_plot <- me1
bottom_plots <- wrap_plots(main_effect_plots, ncol = 2)
combined_main_effects <- top_plot / bottom_plots + plot_annotation(title = "Predicted discrepancy, Path vs AI", theme = theme(plot.title = element_text(size = 18)))

# Save the combined plot
ggsave("combined_main_effects.pdf", combined_main_effects, width = 8, height = 7, dpi = 300)


## fmi



coef_names <- names(classic_no_forced_interactions$TCC_FMI_minus_TCC_AI$model$coefficients)[-1]
# Remove "log_" or "sq_" prefixes and get unique base variables
clean_names <- str_replace(coef_names, "^(log_|sq_)", "")
clean_names <- str_replace(clean_names, "\\.(L|C|Q)$", "")
clean_names <- str_replace(clean_names, "Poor/Satisfactory$", "")
clean_names <- str_replace(clean_names, "High$", "")


unique_clean_names <- unique(clean_names)
# Exclude specific variables
exclude_vars <- c(
    "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cells.in.normal.tissue",
    "Foundation.Model...TME.v1.0.1.alpha...Key.Result...Total.area.of.tumor...mm.2"
)
final_vars <- setdiff(unique_clean_names, exclude_vars)

# Generate plots for each final variable
main_effect_plots_fmi <- lapply(seq_along(final_vars), function(i) {
    var <- final_vars[i]
    ylim <- if (i <= 2) c(-10, 50) else c(-25, 25)
    plot_main_effect(
        model = classic_no_forced_interactions$TCC_FMI_minus_TCC_AI$model,
        model_data = all_model_variables_renamed,
        var_base = var,
        get_pretty_name_v2,
        unscale_center = TRUE,
        y_limits = ylim
    ) + ggtitle(NULL) +
        ylab(NULL)
})

# Optionally, name the list elements for reference
names(main_effect_plots_fmi) <- final_vars

# Save each main effect plot as a PNG file
for (i in seq_along(main_effect_plots_fmi)) {
    var_name <- names(main_effect_plots_fmi)[i]
    file_name <- paste0("fmi_test_", var_name, ".png")
    ggsave(filename = file_name, plot = main_effect_plots_fmi[[i]], width = 8, height = 6, dpi = 300)
}

me3 <- plot_main_effect(
    model = classic_no_forced_interactions$TCC_FMI_minus_TCC_AI$model,
    model_data = all_model_variables_renamed,
    var_base = "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...Total.number.of.cells.in.normal.tissue",
    var_base_2 = "Foundation.Model...TME.v1.0.1.alpha...Key.Result...Total.area.of.tumor...mm.2",
    var_base_2_values = c(1, 10, 100),
    get_pretty_name_v2,
    unscale_center = TRUE
) + theme(legend.position = "bottom") + ggtitle(NULL) +
    ylab(NULL) +
    scale_y_continuous(breaks = seq(-100, 100, 20))

ggsave("test_fmi.png", me3, width = 8, height = 6, dpi = 300)

# Modify the first two plots in main_effect_plots_fmi to add ylim(c(-10, 50))
if (length(main_effect_plots_fmi) >= 1) {
    main_effect_plots_fmi[[1]] <- main_effect_plots_fmi[[1]]
}
if (length(main_effect_plots_fmi) >= 2) {
    main_effect_plots_fmi[[2]] <- main_effect_plots_fmi[[2]]
}

# Combine me3 on top and the four main_effect_plots_fmi in a 2x2 grid underneath
top_plot_fmi <- me3
bottom_plots_fmi <- wrap_plots(main_effect_plots_fmi, ncol = 2)
combined_main_effects_fmi <- top_plot_fmi / bottom_plots_fmi + plot_annotation(title = "Predicted discrepancy, Path vs FMI", theme = theme(plot.title = element_text(size = 18)))

# Save the combined plot
ggsave("combined_main_effects_fmi.pdf", combined_main_effects_fmi, width = 8, height = 9, dpi = 300)
