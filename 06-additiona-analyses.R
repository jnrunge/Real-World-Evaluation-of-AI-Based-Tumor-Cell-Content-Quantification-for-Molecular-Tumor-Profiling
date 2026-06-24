################################################################################
# Script: 06-additiona-analyses.R
# Purpose: Additional analyses - pairwise heatmaps of raw TCC values
# Author: Jan-Niklas Runge
#
# Description:
#   - Uses data_df_complete (contains raw TCC_Patho, TCC_AI, TCC_FMI)
#   - Bins each TCC variable into intervals of 5 percentage points
#   - Creates count heatmaps for all three pairwise combinations
#
# Inputs:
#   - data_df_complete (via 01-load-data.R)
#
# Outputs:
#   - output/additional_analyses/tcc_pairwise_heatmaps.pdf
#
################################################################################

output_dir <- file.path(project_dir, "output/additional_analyses")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

tcc_vars <- c("TCC_Patho", "TCC_FMI", "TCC_AI")

tcc_labels <- c(
  TCC_Patho = "TCC Pathologist (VQ)",
  TCC_FMI   = "TCC FMI (MQ)",
  TCC_AI    = "TCC AI (DQ)"
)

bin_tcc <- function(x, by) {
  breaks <- seq(0, 100, by = by)
  cut(x, breaks = breaks, include.lowest = TRUE, right = FALSE)
}

# Extract the lower bound from a factor level like "[15,20)" or "[95,100]"
bin_lower <- function(lev) {
  as.numeric(
    sub("^[\\[\\(]([0-9.]+),.*", "\\1", as.character(lev), perl = TRUE)
  )
}

# Convert cut() labels like "[15,20)" to "15–20"
prettify_bin <- function(lev) {
  s <- as.character(lev)
  s <- sub("^[\\[\\(]", "", s, perl = TRUE)
  s <- sub("[\\]\\)]$", "", s, perl = TRUE)
  sub(",", "–", s, fixed = TRUE)
}

threshold_cutoff <- 20   # clinical TCC threshold
discrepancy_abs  <- 20   # flag if bin lower bounds differ by >= this many units

pairs <- combn(tcc_vars, 2, simplify = FALSE)

# bin sizes chosen so the 20-unit threshold always falls on a bin boundary
for (bin_size in c(5, 10, 20)) {
  plots <- lapply(pairs, function(pair) {
    x_var <- pair[1]
    y_var <- pair[2]

    df <- data_df_complete %>%
      select(all_of(c(x_var, y_var))) %>%
      rename(x = 1, y = 2) %>%
      filter(!is.na(x), !is.na(y)) %>%
      mutate(
        x_bin = bin_tcc(x, bin_size),
        y_bin = bin_tcc(y, bin_size)
      ) %>%
      count(x_bin, y_bin) %>%
      mutate(
        label = paste0(n, " / ", sprintf("%.1f%%", n / sum(n) * 100)),
        x_lower = bin_lower(x_bin),
        y_lower = bin_lower(y_bin),
        threshold_cross =
          (x_lower < threshold_cutoff) != (y_lower < threshold_cutoff),
        large_discrepancy = abs(x_lower - y_lower) >= discrepancy_abs,
        border_color = case_when(
          threshold_cross   ~ "#e63946",  # red: threshold crossed
          large_discrepancy ~ "#f4a261",  # orange: large discrepancy
          TRUE              ~ "#4caf50"   # green: concordant
        )
      )

    ggplot(df, aes(x = x_bin, y = y_bin)) +
      geom_tile(
        aes(fill = n, color = border_color),
        width = 0.85, height = 0.85, linewidth = 1.0
      ) +
      geom_text(aes(label = label), size = 3.5, color = "white") +
      scale_fill_gradient(low = "#666666", high = "#000000") +
      scale_color_identity() +
      scale_x_discrete(labels = prettify_bin, drop = FALSE) +
      scale_y_discrete(labels = prettify_bin, drop = FALSE) +
      guides(fill = "none") +
      labs(
        x = tcc_labels[x_var],
        y = tcc_labels[y_var],
        caption = paste0(
          "Border: red = one method < 20, other ≥ 20  |  ",
          "orange = bins differ by ≥ ", discrepancy_abs,
          " units  |  green = concordant"
        )
      ) +
      theme_bw(11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.caption = element_text(size = 8, hjust = 0)
      )
  })

  pdf_path <- file.path(
    output_dir,
    paste0("tcc_pairwise_heatmaps_bin", bin_size, ".pdf")
  )
  pdf(pdf_path, width = 9, height = 8)
  for (p in plots) print(p)
  dev.off()

  message("Saved: ", pdf_path)
}

# Clinical binning: 0–9, 10–19, 20–50, >50
clinical_labels <- c("0–9", "10–19", "20–50", ">50")
clinical_lower  <- setNames(c(0, 10, 20, 51), clinical_labels)

bin_clinical <- function(x) {
  cut(
    x,
    breaks = c(0, 10, 20, 51, 101),
    labels = clinical_labels,
    include.lowest = TRUE, right = FALSE
  )
}

plots_clinical <- lapply(pairs, function(pair) {
  x_var <- pair[1]
  y_var <- pair[2]

  df <- data_df_complete %>%
    select(all_of(c(x_var, y_var))) %>%
    rename(x = 1, y = 2) %>%
    filter(!is.na(x), !is.na(y)) %>%
    mutate(
      x_bin = bin_clinical(x),
      y_bin = bin_clinical(y)
    ) %>%
    count(x_bin, y_bin) %>%
    mutate(
      label = paste0(n, " / ", sprintf("%.1f%%", n / sum(n) * 100)),
      x_lower = clinical_lower[as.character(x_bin)],
      y_lower = clinical_lower[as.character(y_bin)],
      threshold_cross =
        (x_lower < threshold_cutoff) != (y_lower < threshold_cutoff),
      large_discrepancy = abs(x_lower - y_lower) >= discrepancy_abs,
      border_color = case_when(
        threshold_cross   ~ "#e63946",
        large_discrepancy ~ "#f4a261",
        TRUE              ~ "#4caf50"
      )
    )

  ggplot(df, aes(x = x_bin, y = y_bin)) +
    geom_tile(
      aes(fill = n, color = border_color),
      width = 0.85, height = 0.85, linewidth = 1.0
    ) +
    geom_text(aes(label = label), size = 3.5, color = "white") +
    scale_fill_gradient(low = "#666666", high = "#000000") +
    scale_color_identity() +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    guides(fill = "none") +
    labs(
      x = tcc_labels[x_var],
      y = tcc_labels[y_var],
      caption = paste0(
        "Border: red = one method < 20, other ≥ 20  |  ",
        "orange = bins differ by ≥ ", discrepancy_abs,
        " units  |  green = concordant"
      )
    ) +
    theme_bw(11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 8, hjust = 0)
    )
})

pdf_path <- file.path(output_dir, "tcc_pairwise_heatmaps_clinical.pdf")
pdf(pdf_path, width = 9, height = 8)
for (p in plots_clinical) print(p)
dev.off()
message("Saved: ", pdf_path)
