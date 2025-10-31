# Plot all columns in all_model_variables: histogram for numeric, barplot for ordered/character
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(forcats)

source("discrepancies/2025-10-Data-Version/00-universal-dependencies.R")

# Helper function to plot a single variable
plot_var <- function(df, var, theme_fun = theme_bw, theme_args = list()) {
    v <- df[[var]]
    p <- NULL
    xlab <- get_pretty_name(var)
    print(var)
    print(xlab)
    print(get_pretty_name)
    if (is.numeric(v)) {
        p <- ggplot(df, aes(x = .data[[var]])) +
            geom_histogram(bins = 30, fill = "#72c2ff", color = "#013d6b", alpha = 0.8) +
            labs(title = NULL, x = xlab, y = "Count")
    } else if (is.ordered(v) || is.factor(v) || is.character(v)) {
        # Convert to factor if not already
        v_fac <- if (is.factor(v)) v else factor(v)
        p <- ggplot(df, aes(x = fct_infreq(v_fac))) +
            geom_bar(fill = "#72c2ff", color = "#013d6b", alpha = 0.8) +
            labs(title = NULL, x = xlab, y = "Count") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    if (!is.null(p)) {
        # Apply theme function and any additional theme arguments
        p <- p + do.call(theme_fun, theme_args)
    }
    p
}

# Choose theme and theme arguments here for all plots
my_theme <- theme_bw
my_theme_args <- list(base_size = 14)

# Generate plots for all columns
all_plots <- lapply(names(data_df_renamed), function(var) {
    plot_var(data_df_renamed, var, theme_fun = my_theme, theme_args = my_theme_args)
})

# Remove NULLs (in case of unsupported types)
all_plots <- Filter(Negate(is.null), all_plots)

# Combine all plots into one grid (patchwork)
# Choose number of columns for the grid
ncol_grid <- 4
combined_plot <- wrap_plots(all_plots, ncol = ncol_grid)

# Show or save the combined plot
combined_plot