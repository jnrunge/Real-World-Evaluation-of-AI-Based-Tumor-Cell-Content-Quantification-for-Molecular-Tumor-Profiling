################################################################################
# Script: 02b-tobi-inspired-plot.R
# Purpose: Tobi-inspired pairplots for TCC comparisons (AI vs Patho vs FMI)
# Author: Jan-Niklas Runge
#
# Description:
#   - Builds pairwise scatter/hist matrices overall and per tissue
#   - Lower panels: scatter 
#   - Upper panels: big colored box with correlation value
#
# Inputs:
#   - Requires environment from 00-universal-dependencies.R
#   - Requires processed data from 01-load-data.R (data_df_complete)
#
# Outputs:
#   - output/tobi_plots/*.PDF and output/tobi_plots/poster/*.PDF
################################################################################

# Configuration ---------------------------------------------------------------
base_size <- 14
base_size_poster <- 22

PANEL_LIMITS <- c(0, 100)
HIST_BINS <- 30
HIST_X_LIMS <- c(-80, 80)
HIST_Y_LIMS <- c(0, 30)

CORR_PALETTE <- c("steelblue", "white", "firebrick")
CORR_DOMAIN   <- c(-1, 1)     # for lower-panel corner box
CORR_BOX_DOMAIN <- c(0.15, 0.85) # for big box in upper panel (more contrast)

# Lower-panel corner box position (in data coordinates of 0..100 axes)
CORNER_XMIN <- 80
CORNER_XMAX <- 98
CORNER_YMIN <- 2
CORNER_YMAX <- 20
CORNER_TEXT_SIZE <- 4.2
CORNER_ALPHA <- 0.75

# Ensure deps and data are available -----------------------------------------
source(file.path(project_dir, "00-universal-dependencies.R"))
if (!exists("data_df_complete")) {
  message("data_df_complete not found. Sourcing 01-load-data.R ...")
  source(file.path(project_dir, "01-load-data.R"))
}

# Helpers ---------------------------------------------------------------------

#' Map correlation to color
get_corr_color <- function(x, y, domain = CORR_DOMAIN) {
  r <- suppressWarnings(cor(x, y, use = "complete.obs"))
  col <- scales::col_numeric(palette = CORR_PALETTE, domain = domain, na.color = "grey80")(r)
  list(r = r, col = col)
}

#' Upper panel: single big colored box with r
add_corr_box_panel <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  if (identical(mapping$x, mapping$y)) {
    return(ggplot() + theme_void())
  }
  r_info <- get_corr_color(x, y, domain = CORR_BOX_DOMAIN)
  sq_xmin <- 15; sq_xmax <- 85; sq_ymin <- 15; sq_ymax <- 85
  sq_size <- sq_xmax - sq_xmin

  ggplot() +
    annotate("rect",
      xmin = sq_xmin, xmax = sq_xmax, ymin = sq_ymin, ymax = sq_ymax,
      fill = r_info$col, alpha = 0.85, color = "black"
    ) +
    annotate("text",
      x = sq_xmin + sq_size/2, y = sq_ymin + sq_size/2,
      label = sprintf("%.2f", r_info$r), size = 8, fontface = "bold"
    ) +
    scale_x_continuous(limits = PANEL_LIMITS) +
    scale_y_continuous(limits = PANEL_LIMITS) +
    theme_void()
}

#' Lower panel: scatter + identity + small correlation corner box
add_identity_line_colored_corner <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  r_info <- get_corr_color(x, y, domain = CORR_DOMAIN)

  ggplot(data = data, mapping = mapping) +
    geom_point(aes(color = `Sample type`), alpha = 0.8, size = 2, shape = 1) +
    geom_abline(intercept = 0, slope = 1, color = "grey60", linetype = "dashed", size = 0.7, alpha = 0.5) +
    # annotate("rect",
    #   xmin = CORNER_XMIN, xmax = CORNER_XMAX, ymin = CORNER_YMIN, ymax = CORNER_YMAX,
    #   fill = r_info$col, alpha = CORNER_ALPHA, color = "black"
    # ) +
    # annotate("text",
    #   x = (CORNER_XMIN + CORNER_XMAX)/2, y = (CORNER_YMIN + CORNER_YMAX)/2,
    #   label = sprintf("r=%.2f", r_info$r), size = CORNER_TEXT_SIZE, fontface = "bold"
    # ) +
    scale_x_continuous(limits = PANEL_LIMITS) +
    scale_y_continuous(limits = PANEL_LIMITS) +
    scale_color_manual(values=c("Biopsy" = "#1B9E77",
                        "Cytology" = "#D95F02",
                        "Resection" = "#7570B3"))
}

#' Apply common scale/theme tweaks to a ggpairs matrix
apply_common_pairplot_tweaks <- function(pmat) {
  for (i in 1:3) {
    for (j in 1:3) {
      if (i != j) {
        pmat[i, j] <- pmat[i, j] +
          scale_x_continuous(limits = PANEL_LIMITS) +
          scale_y_continuous(limits = PANEL_LIMITS) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      } else {
        pmat[i, j] <- pmat[i, j] +
          scale_x_continuous(limits = PANEL_LIMITS) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }
    }
  }
  pmat + theme(
    strip.background = element_rect(fill = "#e3f0fc", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )
}

#' Build a pairplot matrix with custom panels
make_pairplot <- function(df, base_sz, title = NULL) {
  ggpairs(
    df,
    columns = 1:3,
    upper = list(continuous = add_corr_box_panel),
    diag = list(continuous = wrap("barDiag", fill = "#72c2ff", color = "#013d6b", bins = HIST_BINS)),
    lower = list(continuous = add_identity_line_colored_corner),
    columnLabels = c("AI", "Patho", "FMI")
  ) +
    theme_bw(base_sz) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(title %||% "")
}

#' Extract legend as patchwork element
get_legend <- function(p) {
  g <- ggplotGrob(
    p + theme(legend.position = "right") +
      guides(color = guide_legend(override.aes = list(shape = 16, fill = NA)))
  )
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (!length(idx)) stop("No legend found in plot")
  patchwork::wrap_elements(g$grobs[[idx]])
}

# Data selection --------------------------------------------------------------
df <- data_df_complete %>% dplyr::select(TCC_AI, TCC_Patho, TCC_FMI, Cancer, `Sample type`)

# Overall matrix --------------------------------------------------------------
p_matrix <- make_pairplot(df, base_size, title = "All combined")
p_matrix <- apply_common_pairplot_tweaks(p_matrix)

# By-tissue matrices ----------------------------------------------------------
tissue_levels <- unique(df$Cancer)
tissue_levels <- tissue_levels[!is.na(tissue_levels)]

p_matrix_tissue_list <- map(
  tissue_levels,
  function(tissue) {
    df_sub <- df %>% dplyr::filter(Cancer == tissue)
    pm <- make_pairplot(df_sub, base_size, title = paste(tissue, "CA"))
    apply_common_pairplot_tweaks(pm)
  }
)

# By-sample-type matrices
sample_type_levels <- unique(df$`Sample type`)
sample_type_levels <- sample_type_levels[!is.na(sample_type_levels)]

p_matrix_sampletype_list <- map(
  sample_type_levels,
  function(stype) {
    df_sub <- df %>% dplyr::filter(`Sample type` == stype)
    pm <- make_pairplot(df_sub, base_size, title = stype)
    apply_common_pairplot_tweaks(pm)
  }
)

# Combine matrices in grid ----------------------------------------------------
all_pmatrices <- c(list(p_matrix), p_matrix_tissue_list, p_matrix_sampletype_list)
n_plots <- length(all_pmatrices)
ncol_grid <- 3

combined_p_matrix <- wrap_elements(ggmatrix_gtable(all_pmatrices[[1]]))
if (n_plots > 1) {
    for (i in 2:n_plots) {
        combined_p_matrix <- combined_p_matrix | wrap_elements(ggmatrix_gtable(all_pmatrices[[i]]))
        if ((i %% ncol_grid) == 0 && i != n_plots) {
            combined_p_matrix <- combined_p_matrix / NULL
        }
    }
}
if (n_plots == 1) {
    combined_p_matrix <- combined_p_matrix
} else {
    combined_p_matrix <- combined_p_matrix + plot_layout(ncol = ncol_grid, nrow = ceiling(n_plots / ncol_grid))
}

# Legend handling and save ----------------------------------------------------
# Extract a lower panel (with points) from the first p_matrix
point_panel <- p_matrix[2, 1] # or any lower panel with points
legend_sample_type <- get_legend(point_panel)

# Combine: legend on the left, combined_p_matrix on the right
final_combined <- (legend_sample_type | combined_p_matrix) + plot_layout(widths = c(0.08, 0.92))

dir.create(file.path(project_dir, "output/tobi_plots"), showWarnings = FALSE, recursive = TRUE)

# Save the combined plot (optional)
ggsave(file.path(project_dir, "output/tobi_plots/TCC_correlations_by_tissue.PDF"), final_combined, width = 18, height = 10, dpi = 300)

ggsave(file.path(project_dir, "output/tobi_plots/TCC_correlations.PDF"), p_matrix, width = 10, height = 10, dpi = 300)

# Histogram (sample-type discrepancy) ----------------------------------------
sampletype_discrepancy_hist <- ggplot(data_df_complete, aes(x = TCC_Patho_minus_TCC_AI, fill = `Sample type`)) +
  geom_histogram(position = "dodge", bins = HIST_BINS, color = "#013d6b", alpha = 0.85) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 0.5, alpha = 0.5) +
  facet_wrap(~`Sample type`, ncol = 1) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw(base_size) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#e3f0fc", color = "black"),
    axis.text = element_text(size = base_size),
    axis.title = element_text(size = base_size + 2),
    plot.title = element_text(size = base_size + 4, face = "bold"),
    legend.position = "none"
  ) +
  labs(x = "«Path» - «AI»", y = "Frequency (n)", title = NULL) +
  guides(x = guide_axis(minor.ticks = TRUE), y = guide_axis(minor.ticks = TRUE)) +
  scale_x_continuous(limits = HIST_X_LIMS, breaks = round(seq(HIST_X_LIMS[1], HIST_X_LIMS[2], length.out = 5)), minor_breaks = round(seq(HIST_X_LIMS[1], HIST_X_LIMS[2], by = 10))) +
  scale_y_continuous(limits = HIST_Y_LIMS, breaks = seq(HIST_Y_LIMS[1], HIST_Y_LIMS[2], by = 5), minor_breaks = NULL)

ggsave(file.path(project_dir, "output/tobi_plots/TC_Path_minus_TC_AI_by_sample_type.PDF"), sampletype_discrepancy_hist, width = 4, height = 10, dpi = 300)

# Combine final_combined and sampletype_discrepancy_hist into a new plot with sampletype_discrepancy_hist on the right



# Remove legend from sampletype_discrepancy_hist for a cleaner layout
sampletype_discrepancy_hist_noleg <- sampletype_discrepancy_hist + theme(legend.position = "none")

# Combine: final_combined on the left, sampletype_discrepancy_hist on the right
combined_with_hist <- (combined_p_matrix | sampletype_discrepancy_hist_noleg) +
    plot_layout(widths = c(0.85, 0.15))

# Show or save the combined plot
ggsave(file.path(project_dir, "output/tobi_plots/TCC_correlations_with_sampletype_hist.PDF"), combined_with_hist, width = 15, height = 12, dpi = 300)

# Poster variants -------------------------------------------------------------
# Poster overall plot matrix
p_matrix_poster <- make_pairplot(df, base_size_poster, title = "All combined")
p_matrix_poster <- apply_common_pairplot_tweaks(p_matrix_poster)

# Poster faceted-by-tissue
p_matrix_tissue_list_poster <- purrr::map(
  tissue_levels,
  function(tissue) {
    df_sub <- df %>% dplyr::filter(Cancer == tissue)
    pm <- make_pairplot(df_sub, base_size_poster, title = paste(tissue, "CA"))
    apply_common_pairplot_tweaks(pm)
  }
)

# Poster faceted-by-sample-type
sample_type_levels <- unique(df$`Sample type`)
sample_type_levels <- sample_type_levels[!is.na(sample_type_levels)]

p_matrix_sampletype_list_poster <- purrr::map(
  sample_type_levels,
  function(stype) {
    df_sub <- df %>% dplyr::filter(`Sample type` == stype)
    pm <- make_pairplot(df_sub, base_size_poster, title = stype)
    apply_common_pairplot_tweaks(pm)
  }
)

dir.create(file.path(project_dir, "output/tobi_plots/poster"), showWarnings = FALSE, recursive = TRUE)

# Poster combined matrices
all_pmatrices_poster <- c(list(p_matrix_poster), p_matrix_tissue_list_poster, p_matrix_sampletype_list_poster)
n_plots_poster <- length(all_pmatrices_poster)
ncol_grid_poster <- 3

combined_p_matrix_poster <- patchwork::wrap_elements(GGally::ggmatrix_gtable(all_pmatrices_poster[[1]]))
if (n_plots_poster > 1) {
    for (i in 2:n_plots_poster) {
        combined_p_matrix_poster <- combined_p_matrix_poster | patchwork::wrap_elements(GGally::ggmatrix_gtable(all_pmatrices_poster[[i]]))
        if ((i %% ncol_grid_poster) == 0 && i != n_plots_poster) {
            combined_p_matrix_poster <- combined_p_matrix_poster / NULL
        }
    }
}
if (n_plots_poster == 1) {
    combined_p_matrix_poster <- combined_p_matrix_poster
} else {
    combined_p_matrix_poster <- combined_p_matrix_poster + patchwork::plot_layout(ncol = ncol_grid_poster, nrow = ceiling(n_plots_poster / ncol_grid_poster))
}

# Poster histogram
sampletype_discrepancy_hist_poster <- ggplot(data_df_complete, aes(x = TCC_Patho_minus_TCC_AI, fill = `Sample type`)) +
    geom_histogram(position = "dodge", bins = HIST_BINS, color = "#013d6b", alpha = 0.85) +
    facet_wrap(~`Sample type`, ncol = 1) +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw(base_size_poster) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#e3f0fc", color = "black"),
        axis.text = element_text(size = base_size_poster),
        axis.title = element_text(size = base_size_poster + 2),
        plot.title = element_text(size = base_size_poster + 4, face = "bold"),
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
    scale_x_continuous(limits = HIST_X_LIMS, breaks = round(seq(HIST_X_LIMS[1], HIST_X_LIMS[2], length.out = 5)), minor_breaks = round(seq(HIST_X_LIMS[1], HIST_X_LIMS[2], by = 10))) +
    scale_y_continuous(limits = HIST_Y_LIMS, breaks = seq(HIST_Y_LIMS[1], HIST_Y_LIMS[2], by = 5), minor_breaks = NULL)

# Poster combined with histogram (no legend)
combined_with_hist_poster <- (combined_p_matrix_poster | sampletype_discrepancy_hist_poster) +
    patchwork::plot_layout(widths = c(0.85, 0.15))

ggsave(file.path(project_dir, "output/tobi_plots/poster/TCC_correlations_with_sampletype_hist.PDF"),
       combined_with_hist_poster, width = 18, height = 15, dpi = 300)

