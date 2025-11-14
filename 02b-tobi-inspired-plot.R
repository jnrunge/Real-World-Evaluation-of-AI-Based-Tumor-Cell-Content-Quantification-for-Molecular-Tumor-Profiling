source(file.path(project_dir, "00-universal-dependencies.R"))

# pretty & simple tobi-inspired plots
# Pairwise scatterplots and histograms for TC_AI, TC_Path, TC_FMI

# Select relevant columns
df <- data_df_complete %>% dplyr::select(TCC_AI, TCC_Patho, TCC_FMI, Cancer, `Sample type`)

# Function to compute correlation and map to color
get_corr_color <- function(x, y) {
    r <- cor(x, y, use = "complete.obs")
    # Map correlation to color: blue (0), white (0.5), red (1)
    col <- scales::col_numeric(
        palette = c("steelblue", "white", "firebrick"),
        domain = c(0, 1)
    )(r)
    list(r = r, col = col)
}

# Custom plot function for upper/lower panels with squared colored box in bottom right
add_identity_line_colored_corner <- function(data, mapping, ...) {
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    sample_type <- data$`Sample type`
    corr_info <- get_corr_color(x, y)
    # Square size and position

    ggplot(data = data, mapping = mapping) +
        geom_point(aes(color = `Sample type`), alpha = 0.8, size = 2, shape = 1) +
        geom_abline(intercept = 0, slope = 1, color = "grey60", linetype = "dashed", size = 0.7, alpha = 0.5) +
        # Add colored square in bottom right
        # Add correlation text inside the square
        scale_x_continuous(limits = c(0, 100)) +
        scale_y_continuous(limits = c(0, 100)) +
        scale_color_brewer(palette = "Dark2")
}



# Function to compute correlation and map to color
get_r_color <- function(x, y) {
    r <- cor(x, y, use = "complete.obs")
    # Map R² to color: blue (0), white (0.5), red (1)
    col <- scales::col_numeric(
        palette = c("steelblue", "white", "firebrick"),
        domain = c(0.15, 0.85)
    )(r)
    list(r = r, col = col)
}

# Custom upper panel: only show R² box in top right
add_r2_box_panel <- function(data, mapping, ...) {
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    if (identical(mapping$x, mapping$y)) {
        return(ggplot() +
            theme_void())
    }
    r_info <- get_r_color(x, y)
    # Box size and position

    sq_xmin <- 15
    sq_xmax <- 85
    sq_ymin <- 15
    sq_ymax <- 85

    sq_size <- sq_xmax - sq_xmin

    ggplot() +
        annotate("rect",
            xmin = sq_xmin, xmax = sq_xmax, ymin = sq_ymin, ymax = sq_ymax,
            fill = r_info$col, alpha = 0.85, color = "black"
        ) +
        annotate("text",
            x = sq_xmax - sq_size / 2, y = sq_ymax - sq_size / 2,
            label = sprintf("%.2f", r_info$r), size = 8, fontface = "bold"
        ) +
        scale_x_continuous(limits = c(0, 100)) +
        scale_y_continuous(limits = c(0, 100)) +
        theme_void()
}

# 1. The original plot (all tissues together)
p_matrix <- ggpairs(
    df,
    columns = 1:3,
    upper = list(continuous = add_r2_box_panel),
    diag = list(continuous = wrap("barDiag", fill = "#72c2ff", color = "#013d6b", bins = 30)),
    lower = list(continuous = add_identity_line_colored_corner),
    columnLabels = c("AI", "Patho", "FMI")
) +
    theme_bw(14) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle("All combined")

# ...rest of your code unchanged...

for (i in 1:3) {
    for (j in 1:3) {
        if (i != j) {
            p_matrix[i, j] <- p_matrix[i, j] +
                scale_x_continuous(limits = c(0, 100)) +
                scale_y_continuous(limits = c(0, 100)) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        } else {
            p_matrix[i, j] <- p_matrix[i, j] +
                scale_x_continuous(limits = c(0, 100)) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }
    }
}
p_matrix <- p_matrix + theme(strip.background = element_rect(fill = "#e3f0fc", color = "black"))
p_matrix <- p_matrix + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
)

# 2. Faceted by Tissue (split by Tissue)
tissue_levels <- unique(df$Cancer)
tissue_levels <- tissue_levels[!is.na(tissue_levels)]
p_matrix_tissue_list <- map(
    tissue_levels,
    function(tissue) {
        df_sub <- df %>% filter(Cancer == tissue)
        pm <- ggpairs(
            df_sub,
            columns = 1:3,
            upper = list(continuous = add_r2_box_panel),
            diag = list(continuous = wrap("barDiag", fill = "#72c2ff", color = "#013d6b", bins = 30)),
            lower = list(continuous = add_identity_line_colored_corner),
            columnLabels = c("AI", "Patho", "FMI")
        ) +
            theme_bw(14) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            ggtitle(paste(tissue, "CA"))
        for (i in 1:3) {
            for (j in 1:3) {
                if (i != j) {
                    pm[i, j] <- pm[i, j] +
                        scale_x_continuous(limits = c(0, 100)) +
                        scale_y_continuous(limits = c(0, 100)) +
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                } else {
                    pm[i, j] <- pm[i, j] +
                        scale_x_continuous(limits = c(0, 100)) +
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                }
            }
        }
        pm <- pm + theme(strip.background = element_rect(fill = "#e3f0fc", color = "black"))
        pm <- pm + theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank()
        )
        pm
    }
)

# Optionally, combine all plots (original + by tissue) using patchwork
# Combine all plots (original + by tissue) in a grid with ncol = ceiling((1 + length(p_matrix_tissue_list)) / 2)
all_pmatrices <- c(list(p_matrix), p_matrix_tissue_list)
n_plots <- length(all_pmatrices)
ncol_grid <- ceiling(n_plots / 2)

combined_p_matrix <- wrap_elements(ggmatrix_gtable(all_pmatrices[[1]]))
if (n_plots > 1) {
    for (i in 2:n_plots) {
        combined_p_matrix <- combined_p_matrix | wrap_elements(ggmatrix_gtable(all_pmatrices[[i]]))
        if ((i %% ncol_grid) == 0 && i != n_plots) {
            combined_p_matrix <- combined_p_matrix / NULL
        }
    }
}
# If only one row, ensure it's a patchwork object
if (n_plots == 1) {
    combined_p_matrix <- combined_p_matrix
} else {
    combined_p_matrix <- combined_p_matrix + plot_layout(ncol = ncol_grid, nrow = ceiling(n_plots / ncol_grid))
}

get_legend <- function(p) {
    # Force legend to use filled dots (shape = 16) for the legend only
    g <- ggplotGrob(
        p + theme(legend.position = "right") +
            guides(color = guide_legend(override.aes = list(shape = 16, fill = NA)))
    )
    legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
    if (length(legend_index) == 0) stop("No legend found in plot")
    patchwork::wrap_elements(g$grobs[[legend_index]])
}

# Extract a lower panel (with points) from the first p_matrix
point_panel <- p_matrix[2, 1] # or any lower panel with points
legend_sample_type <- get_legend(point_panel)

# Combine: legend on the left, combined_p_matrix on the right
final_combined <- (legend_sample_type | combined_p_matrix) + plot_layout(widths = c(0.08, 0.92))

dir.create(file.path(project_dir, "output/tobi_plots"), showWarnings = FALSE)

# Save the combined plot (optional)
ggsave(file.path(project_dir, "output/tobi_plots/TCC_correlations_by_tissue.PDF"), final_combined, width = 18, height = 10, dpi = 300)

ggsave(file.path(project_dir, "output/tobi_plots/TCC_correlations.PDF"), p_matrix, width = 10, height = 10, dpi = 300)

sampletype_discrepancy_hist <- ggplot(data_df_complete, aes(x = TCC_Patho_minus_TCC_AI, fill = `Sample type`)) +
    geom_histogram(position = "dodge", bins = 30, color = "#013d6b", alpha = 0.85) +
    facet_wrap(~`Sample type`, ncol = 1) +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw(14) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#e3f0fc", color = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold"),
        legend.position = "none"
    ) +
    labs(
        x = "«Path» - «AI»",
        y = "Frequency (n)",
        title = NULL
    ) +
    guides(
        x = guide_axis(minor.ticks = TRUE),
        y = guide_axis(minor.ticks = TRUE)
    ) +
    scale_x_continuous(limits = c(-80, 80), breaks = round(seq(-80, 80, length.out = 5)), minor_breaks = round(seq(-80, 80, by = 10))) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5), minor_breaks = NULL)

ggsave(file.path(project_dir, "output/tobi_plots/TC_Path_minus_TC_AI_by_sample_type.PDF"), sampletype_discrepancy_hist, width = 4, height = 10, dpi = 300)

# Combine final_combined and sampletype_discrepancy_hist into a new plot with sampletype_discrepancy_hist on the right



# Remove legend from sampletype_discrepancy_hist for a cleaner layout
sampletype_discrepancy_hist_noleg <- sampletype_discrepancy_hist + theme(legend.position = "none")

# Combine: final_combined on the left, sampletype_discrepancy_hist on the right
combined_with_hist <- (combined_p_matrix | sampletype_discrepancy_hist_noleg) +
    plot_layout(widths = c(0.85, 0.15))

# Show or save the combined plot
ggsave(file.path(project_dir, "output/tobi_plots/TCC_correlations_with_sampletype_hist.PDF"), combined_with_hist, width = 15, height = 9, dpi = 300)

