################################################################################
# Script: 04a-plot-models.R
# Purpose: Generate visualizations for best linear models
# Author: Jan-Niklas Runge
# 
# Description:
#   - Creates coefficient plots
#   - Generates main effect and interaction plots for all model terms
#   - Handles numeric, factor, and interaction terms appropriately
# 
# Outputs:
#   - output/model_plots/TCC_discrepancy_model_*.pdf: Coefficient plots with R²
#   - output/model_plots/*_main_effect_*.pdf: Main effect plots
#   - output/model_plots/*_interaction_*.pdf: Interaction plots
#   - output/model_plots/poster/*.pdf: High-resolution poster versions
# 
# Dependencies: See 00-universal-dependencies.R for package requirements
################################################################################

# Configuration ----------------------------------------------------------------
# Plot font sizes
base_size <- 14
base_size_poster <- 22

# Plot dimensions (base values; adjusted per plot type)
DEFAULT_WIDTH <- 3
DEFAULT_HEIGHT <- 3
COEFF_PLOT_WIDTH <- 10
COEFF_PLOT_HEIGHT <- 13
INTERACTION_PLOT_WIDTH <- 11
INTERACTION_PLOT_HEIGHT <- 8
PLOT_DPI <- 300

# Y-axis limits for discrepancy predictions
DISCREPANCY_YLIM <- c(-100, 100)

# Number of grouping levels for continuous variables in interaction plots
N_GROUPS <- 4

# Color palettes
# Color-blind friendly palette (Okabe-Ito)
CB_PALETTE <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Annotator colors: Intercept (grey), AI (blue), Pathologist (purple)
SCALE_INTERCEPT_AI_PATH <- c("Intercept" = "#797979", "FALSE" = "#A767FF", "TRUE" = "#1041FF")

# Create output directories ----------------------------------------------------
output_dir <- file.path(project_dir, "output/model_plots/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_dir_poster <- file.path(project_dir, "output/model_plots/poster/")
dir.create(output_dir_poster, showWarnings = FALSE, recursive = TRUE)

# Helper functions -------------------------------------------------------------

#' Plot main effect or interaction for a model term
#' 
#' @param model Linear model object
#' @param model_data Scaled model data
#' @param var_base Character string, primary variable name
#' @param get_pretty_name_v2 Function to convert variable names
#' @param ylab Y-axis label (default: "Predicted Discrepancy")
#' @param title_prefix Title prefix (default: "Main Effect of")
#' @param var_base_2 Character string, optional grouping variable
#' @param var_base_2_values Numeric/character vector, values for grouping variable
#' @param palette Color palette for grouping variable
#' @param var_base_2_correction Character string, optional correction variable
#' @param var_base_2_correction_values Numeric vector, correction values
#' @param unscale_center Logical, use raw (unscaled) data for x-axis
#' @param y_limits Numeric vector of length 2, y-axis limits
#' @return ggplot object
#' @details 
#'   - Numeric var_base: line plot with confidence ribbon
#'   - Factor var_base: boxplot of predictive distribution
#'   - Interaction: multiple lines/boxplots colored by var_base_2
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
    var_is_factor <- is.factor(model_data[[var_base]]) || is.character(model_data[[var_base]])
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
            # ensure that val2 is only "0" for level that is always "0"
            for(lvl in levels_base){
                if(var(model_data[[var_base_2]][model_data[[var_base]]==lvl])==0){
                    combos$val2[combos$level==lvl] <- unique(model_data[[var_base_2]][model_data[[var_base]]==lvl])
                    if(unscale_center){
                        combos$val2[combos$level==lvl] <- unique(raw_data[[var_base_2]][raw_data[[var_base]]==lvl])
                    }
            }
            }
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
                } else if (paste0("sq_", var_base_2) %in% model_vars) {
                    nd[[paste0("sq_", var_base_2)]] <- if (unscale_center) {
                        apply_scaler(var2_val^2, s_sq_base2)
                    } else {
                        var2_val^2
                    }
                } else {
                    if (unscale_center) {
                        apply_scaler(var2_val, s_base2)
                    } else {
                        var2_val
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
                x = {
                    lbl <- get_pretty_name_v2(var_base)
                    spaces <- gregexpr(" ", lbl, fixed = TRUE)[[1]]
                    if (spaces[1] != -1) {
                    mid <- nchar(lbl) / 2
                    space_idx <- spaces[which.min(abs(spaces - mid))]
                    substr(lbl, space_idx, space_idx) <- "\n"
                    }
                    lbl
                },
                y = ylab,
                title = paste(
                    title_prefix,
                    get_pretty_name_v2(var_base),
                    if (!is.null(var_base_2)) paste("by", get_pretty_name_v2(var_base_2)) else ""
                )
            ) +
            ggplot2::theme_bw(base_size) +
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
                x = {
                    lbl <- get_pretty_name_v2(var_base)
                    spaces <- gregexpr(" ", lbl, fixed = TRUE)[[1]]
                    if (spaces[1] != -1) {
                    mid <- nchar(lbl) / 2
                    space_idx <- spaces[which.min(abs(spaces - mid))]
                    substr(lbl, space_idx, space_idx) <- "\n"
                    }
                    lbl
                },
                y = ylab,
                title = paste(title_prefix, get_pretty_name_v2(var_base))
            ) +
            scale_x_continuous(labels = scales::comma) +
            theme_bw(base_size) +
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
        theme_bw(base_size) +
        coord_cartesian(ylim = c(-100, 100))
    return(p)
}

#' Convert variable names to pretty display names (extended version)
#' 
#' @param var Character string, variable name
#' @return Character string, formatted display name
#' @details 
#'   - Handles .L/.C/.Q suffixes for ordered factors (linear/cubic/quadratic contrasts)
#'   - Processes interaction terms
#'   - Uses pretty_names lookup table
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

#' Plot model coefficients with pretty names, AI/Path coloring, and R² annotation
#' 
#' @param model Linear model object
#' @param importance_summary Importance summary from best_models output
#' @param get_pretty_name_v2 Function to convert variable names
#' @param title Character string, plot title
#' @return ggplot object
#' @details
#'   - Orders coefficients by importance
#'   - Colors by source: Intercept (grey), AI (blue), Pathologist (purple)
#'   - Adds R² box 
#'   - Symmetric y-axis limits around zero
plot_model_coefficients <- function(model, importance_summary, get_pretty_name_v2, title = "TCC Discrepancy Model") {
    coef_names <- names(model$coefficients[-1])
    orig_vars <- importance_summary$orig_vars

    # Remove .L/.Q/.C suffixes for matching
    coef_base <- sub("\\.(L|C|Q)$", "", coef_names)

    # Remove factor level suffixes like "High" or "Good" from coef_base and orig_vars for matching
    remove_level_suffix <- function(x) {
        suffixes <- c("High", "Good", "Cytology", "Biopsy", "Resection", "Moderate", "Low", "Poor", "Excellent", "Fair")
        pattern <- paste0("(", paste(suffixes, collapse = "|"), ")$")
        process_one <- function(elem) {
            parts <- strsplit(elem, ":", fixed = TRUE)[[1]]
            parts <- vapply(parts, function(p) sub(pattern, "", p), character(1))
            paste(parts, collapse = ":")
        }
        vapply(x, process_one, character(1), USE.NAMES = FALSE)
    }
    coef_base <- remove_level_suffix(coef_base)



    # For each orig_var, find matching coef indices (with or without .L/.Q)
    ordered_indices <- unlist(lapply(orig_vars, function(var) which(coef_base == var)))
    
    # Handle interaction terms: place them right after the lower-ranking component
    interaction_indices <- which(grepl(":", coef_base))
    if (length(interaction_indices) > 0) {
        for (int_idx in interaction_indices) {
            # Get the interaction term
            int_term <- coef_base[int_idx]
            # Split into components
            components <- strsplit(int_term, ":", fixed = TRUE)[[1]]
            
            # Find the positions of each component in ordered_indices
            comp_positions <- sapply(components, function(comp) {
                which(ordered_indices==(which(coef_base==comp) %>% max()))
            })
            comp_positions <- comp_positions[!is.na(comp_positions)]
            
            if (length(comp_positions) > 0) {
                # Place interaction right after the lower-ranking (higher position number) component
                insert_after <- max(comp_positions)
                # Insert the interaction index at this position
                ordered_indices <- append(ordered_indices, int_idx, after = insert_after)
            } else {
                # If components not found in ordered_indices, append at end
                ordered_indices <- c(ordered_indices, int_idx)
            }
        }
    }
    ordered_indices <- ordered_indices[ordered_indices > 0]

    # Compute y-axis limits so that min and max are symmetric around zero
    coef_vals <- model$coefficients
    ordered_indices <- ordered_indices[ordered_indices %in% (which(!is.na(coef_vals))-1)]
    coef_names <- coef_names[which(!is.na(coef_vals))-1] # -1 for intercept already removed
    coef_vals <- coef_vals[!is.na(coef_vals)]
    # Get standard errors for coefficients
    coef_se <- summary(model)$coefficients[, "Std. Error"]
    # Compute range including SE
    coef_range <- range(coef_vals + 1.96 * coef_se, coef_vals - 1.96 * coef_se, na.rm = TRUE)
    abs_max <- max(abs(coef_range))
    y_limits <- c(-abs_max - 5, abs_max + 5)


    # after all that messing around, make sure we dont skip numbers here
    while (length(ordered_indices) < (max(ordered_indices) - min(ordered_indices) + 1)) {
        current_max <- max(ordered_indices)
        if (!(current_max - 1) %in% ordered_indices) {
            ordered_indices[ordered_indices == current_max] <- current_max - 1
        } else {
            break # Safety break if we can't compress further
        }
    }


    pm <- plot_model(
        model,
        show.intercept = TRUE, order.terms = c(1, (ordered_indices + 1)),
        axis.labels = c("(Intercept)" = "Intercept", vapply(
            coef_names,
            get_pretty_name_v2,
            character(1)
        ))
    ) + theme_bw(base_size) + ylab("Coefficient Estimate") + ggtitle(title) +
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


    pm <- pm + theme_bw(base_size) + 
    theme(legend.position = "bottom")
    pm
}

#' Save both regular and poster versions of a plot
#' 
#' @param plot_obj ggplot object
#' @param filename Character string, output filename
#' @param width Numeric, plot width in inches
#' @param height Numeric, plot height in inches
#' @param width_poster Numeric, poster width (defaults to width)
#' @param height_poster Numeric, poster height (defaults to height)
#' @param dpi Numeric, resolution
save_plot_both <- function(plot_obj, filename, width, height, width_poster = NULL, height_poster = NULL, dpi = 300) {
    # Use regular dimensions for poster if not specified
    if (is.null(width_poster)) width_poster <- width
    if (is.null(height_poster)) height_poster <- height
    
    # Regular version
    ggsave(paste0(output_dir, filename), plot_obj, width = width, height = height, dpi = dpi)
    
    # Poster version - increase base_size
    plot_poster <- plot_obj + theme_bw(base_size_poster)+ theme(legend.position = "bottom")
    ggsave(paste0(output_dir_poster, filename), plot_poster, width = width_poster, height = height_poster, dpi = dpi)
}

#' Get sensible grouping values for a variable
#' 
#' @param var_name Character string, variable name
#' @param data Data frame containing variable
#' @param n_groups Integer, number of groups to create
#' @return Character or numeric vector of grouping values
#' @details
#'   - Factors: returns all levels or evenly sampled subset
#'   - Numeric: returns quantiles (10th to 90th percentile)
#'   - Rounds numeric values to sensible precision
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
        if(length(unique(quantiles)) < n_groups){
            if(max(var_data, na.rm=TRUE) == 1 && min(var_data, na.rm=TRUE) == 0){
                # Special case for proportions between 0 and 1
                quantiles <- seq(0, 1, length.out = n_groups)
            } else {
                # Fallback to equally spaced values across range
                quantiles <- seq(min(var_data, na.rm=TRUE), max(var_data, na.rm=TRUE), length.out = n_groups)
            }
        }
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

# Generate coefficient plots ---------------------------------------------------
message("\n=== Generating coefficient plots ===")

# Path vs AI model
pm <- plot_model_coefficients(
    classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$model,
    classic_no_forced_interactions$TCC_Patho_minus_TCC_AI$importance_summary,
    get_pretty_name_v2,
    title = "TCC Discrepancy Pathologist vs AI"
)
save_plot_both(pm, "TCC_discrepancy_model_path_vs_ai.pdf", width = COEFF_PLOT_WIDTH, height = COEFF_PLOT_HEIGHT, width_poster=10, height_poster=12)

# FMI vs AI model
pm_fmi_vs_ai <- plot_model_coefficients(
    classic_no_forced_interactions$TCC_FMI_minus_TCC_AI$model,
    classic_no_forced_interactions$TCC_FMI_minus_TCC_AI$importance_summary,
    get_pretty_name_v2,
    title = "TCC Discrepancy FMI vs AI"
)
save_plot_both(pm_fmi_vs_ai, "TCC_discrepancy_model_fmi_vs_ai.pdf", width = 12, height = 16, width_poster=10, height_poster=21)

# Path vs FMI model
pm_path_vs_fmi <- plot_model_coefficients(
    classic_no_forced_interactions$TCC_Patho_minus_TCC_FMI$model,
    classic_no_forced_interactions$TCC_Patho_minus_TCC_FMI$importance_summary,
    get_pretty_name_v2,
    title = "TCC Discrepancy Pathologist vs FMI"
)
save_plot_both(pm_path_vs_fmi, "TCC_discrepancy_model_path_vs_fmi.pdf", width = COEFF_PLOT_WIDTH, height = COEFF_PLOT_HEIGHT, width_poster=10, height_poster=12)

# Generate main effect and interaction plots ----------------------------------
message("\n=== Generating effect plots for all models ===")

# Loop through all models in classic_no_forced_interactions
for (model_name in names(classic_no_forced_interactions)) {
    message(paste("\nProcessing model:", model_name))
    
    model_obj <- classic_no_forced_interactions[[model_name]]$model
    model_data <- data_df_renamed
    
    # Get predictor names (excluding intercept)
    predictor_names <- names(model_obj$coefficients)[-1]
    
    # Helper: extract base variable name (remove suffixes)
    # Rationale: Handle polynomial contrasts (.L/.Q/.C) and factor levels
    get_base_var <- function(pred_name) {
        # Split on ":" for interaction terms
        parts <- strsplit(pred_name, ":", fixed = TRUE)[[1]]
        
        # Process each part
        parts_clean <- vapply(parts, function(part) {
            # Remove polynomial contrast suffixes
            part_clean <- sub("\\.(L|Q|C)$", "", part)
            # Remove common factor level suffixes
            part_clean <- sub("(High|Good|Low|Bad|Medium|Cytology|Resection|Biopsy)$", "", part_clean)
            return(part_clean)
        }, character(1), USE.NAMES = FALSE)
        
        # Paste back together with ":"
        paste(parts_clean, collapse = ":")
    }
    
    # Track processed variables to avoid duplicates
    processed_vars <- character(0)
    
    for (pred in predictor_names) {
        # Skip if already processed (handles .L, .Q, .C variants)
        base_pred <- get_base_var(pred)
        if (base_pred %in% processed_vars) next
        
        # Check if interaction term
        if (grepl(":", pred, fixed = TRUE)) {
            # Interaction term --------------------------------------------------
            message(paste("  Plotting interaction:", base_pred))
            
            parts <- strsplit(base_pred, ":", fixed = TRUE)[[1]]
            if (length(parts) != 2) next
            
            var_base <- parts[1]
            var_base_2 <- parts[2]
            
            # Get appropriate grouping values for var_base_2
            raw_data <- data_df_pre_scaling
            names(raw_data) <- make.names(names(raw_data), unique = TRUE)
            var_base_2_values <- get_grouping_values(var_base_2, raw_data, n_groups = N_GROUPS)
            
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
                    ggtitle(NULL) +
                    ylab(NULL)
                
                # Create safe filename
                safe_filename <- paste0(
                    model_name, "_interaction_",
                    make.names(var_base), "_by_",
                    make.names(var_base_2), ".pdf"
                )
                
                save_plot_both(p, safe_filename, width = INTERACTION_PLOT_WIDTH, height = INTERACTION_PLOT_HEIGHT, width_poster=12)
                message(paste("    Saved:", safe_filename))
            }, error = function(e) {
                message(paste("    Error plotting interaction", var_base, ":", var_base_2, "-", e$message))
            })
            
            processed_vars <- c(processed_vars, base_pred)
            
        } else {
            # Main effect term --------------------------------------------------
            message(paste("  Plotting main effect:", base_pred))
            
            var_base <- base_pred
            
            tryCatch({
                p <- plot_main_effect(
                    model = model_obj,
                    model_data = model_data,
                    var_base = var_base,
                    get_pretty_name_v2 = get_pretty_name_v2,
                    unscale_center = TRUE
                ) + 
                    ggtitle(NULL) +
                    ylab(NULL)
                
                # Create safe filename
                safe_filename <- paste0(
                    model_name, "_main_effect_",
                    make.names(var_base), ".pdf"
                )
                
                save_plot_both(p, safe_filename, width = DEFAULT_WIDTH, height = DEFAULT_HEIGHT, width_poster=12)
                message(paste("    Saved:", safe_filename))
            }, error = function(e) {
                message(paste("    Error plotting main effect", var_base, "-", e$message))
            })
            
            processed_vars <- c(processed_vars, var_base)
        }
    }
}

message("\n=== Model plotting complete ===")
message(paste("Plots saved to:", output_dir))
message(paste("Poster plots saved to:", output_dir_poster))
