library(vip)
library(sjPlot)
pretty_names <- read_csv("discrepancies/variable_pretty_names.csv")
source("discrepancies/2025-10-Data-Version/functions/best_models.R")
source("discrepancies/2025-10-Data-Version/00-universal-dependencies.R")
response_vars <- variables %>%
    filter(type == "response") %>%
    pull(variable)

rds_file <- "discrepancies/2025-10-Data-Version/processed_data/best_models_classic_no_forced_interactions.rds"

if (file.exists(rds_file)) {
    classic_no_forced_interactions <- readRDS(rds_file)
} else {
    classic_no_forced_interactions <- list()
    
    for (response_var in response_vars) {
        var_table <- variables %>%
            filter(type != "response" | variable == response_var)
        
        classic_no_forced_interactions[[response_var]] <- perform_forward_selection(var_table %>% filter(type == "response" | !grepl("TCC", variable)), data_df, force_interaction = FALSE, steps = 5000, seed = 1, reps = 32)
    }
    
    saveRDS(classic_no_forced_interactions, file = rds_file)

    # descrribe alternative models
    # Extract predictors from models and create presence/absence table
    predictor_presence <- list()

    for (response_var in response_vars) {
        # Extract all predictors from the 32 fits
        all_predictors <- lapply(1:32, function(i) {
            model <- classic_no_forced_interactions[[response_var]]$fits[[i]]$model
            # Get predictor names (excluding intercept and response variable)
            pred_names <- names(coef(model))[-1]
            pred_names
        })

        chosen_fit <- classic_no_forced_interactions[[response_var]]$best_idx

        # Get unique predictors across all fits
        unique_predictors <- unique(unlist(all_predictors))

        # Create presence/absence matrix
        presence_matrix <- matrix(0, nrow = length(unique_predictors), ncol = 32)
        rownames(presence_matrix) <- unique_predictors
        colnames(presence_matrix) <- paste0("Fit_", 1:32)

        for (i in 1:32) {
            presence_matrix[all_predictors[[i]], i] <- 1
        }

        predictor_presence[[response_var]] <- as.data.frame(presence_matrix)

        # Perform hierarchical clustering on fits
        # Calculate distance between fits based on predictor presence
        fit_dist <- dist(t(presence_matrix), method = "binary")
        fit_hclust <- hclust(fit_dist, method = "complete")

        # Get dendrogram order
        fit_order <- fit_hclust$order

        # Identify groups of identical fits
        fit_groups <- character(32)
        group_id <- 1
        for (i in 1:32) {
            if (fit_groups[i] == "") {
                # This fit hasn't been assigned to a group yet
                fit_groups[i] <- as.character(group_id)
                # Find all other fits with identical predictor sets
                for (j in (i + 1):32) {
                    if (j <= 32 && all(presence_matrix[, i] == presence_matrix[, j])) {
                        fit_groups[j] <- as.character(group_id)
                    }
                }
                group_id <- group_id + 1
            }
        }

        # Create a data frame for plotting with group information
        plot_data <- predictor_presence[[response_var]] %>%
            mutate(Predictor = rownames(.)) %>%
            pivot_longer(cols = starts_with("Fit_"), names_to = "Fit", values_to = "Present") %>%
            mutate(
                Fit_num = as.numeric(gsub("Fit_", "", Fit)),
                Fit = factor(Fit, levels = paste0("Fit_", fit_order)),
                is_chosen = (Fit == paste0("Fit_", chosen_fit)),
                fit_group = fit_groups[Fit_num]
            )

        # Create dendrogram plot
        library(ggdendro)
        dend_data <- dendro_data(fit_hclust)
        dend_plot <- ggplot(segment(dend_data)) +
            geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
            theme_void() +
            theme(plot.margin = margin(5, 5, 0, 5))


        plot_data$Predictor <- sapply(plot_data$Predictor, function(var) {
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
        }, USE.NAMES = FALSE) %>%
            lapply(get_pretty_name) %>%
            unlist()

        # Plot presence/absence
        main_plot <- plot_data %>%
            ggplot(aes(x = Fit, y = Predictor, fill = interaction(factor(Present), fit_group))) +
            geom_tile(aes(color = is_chosen, linewidth = is_chosen)) +
            scale_fill_manual(values = c(
                setNames(
                    rep("white", length(unique(fit_groups))),
                    paste0("0.", unique(fit_groups))
                ),
                setNames(
                    rainbow(length(unique(fit_groups)), s = 0.7, v = 0.8),
                    paste0("1.", unique(fit_groups))
                )
            ), guide = "none") +
            scale_color_manual(values = c("FALSE" = "white", "TRUE" = "red"), guide = "none") +
            scale_linewidth_manual(values = c("FALSE" = 0.5, "TRUE" = 2), guide = "none") +
            labs(
                title = paste("Predictor Presence Across Model Fits:", response_var),
                subtitle = paste(
                    "Red border indicates chosen fit:", chosen_fit,
                    "| Identical fits have same color |",
                    length(unique(fit_groups)), "unique model(s)"
                ),
                x = "Model Fit (ordered by similarity)", y = "Predictor", fill = ""
            ) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

        # Combine dendrogram and main plot
        library(patchwork)
        combined_plot <- dend_plot / main_plot +
            plot_layout(heights = c(1, 4))

        ggsave(
            paste0(
                "discrepancies/2025-10-Data-Version/model_plots/predictor_presence_",
                response_var, ".pdf"
            ),
            plot = combined_plot,
            width = 12, height = 10
        )
    }
}





## interpret models 
# use compare_model_drop_terms(model, terms) to get term importance/significance