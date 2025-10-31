n_boot <- 1000

library(foreach)
library(doParallel)
library(dplyr)
source("discrepancies/manual-stepwise.R")
source("discrepancies/2025-10-Data-Version/functions/bootstrapping_models.R")
response_vars <- variables %>%
    filter(type == "response") %>%
    pull(variable)

n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

predictor_vars <- variables %>%
    filter(type != "response") %>%
    pull(variable)

available_predictors <- intersect(predictor_vars, names(data_df))

corr_data <- data_df %>%
    dplyr::select(any_of(available_predictors))

ordered_cols <- names(corr_data)[vapply(corr_data, is.ordered, logical(1))]
numeric_cols <- names(corr_data)[vapply(corr_data, is.numeric, logical(1))]
corr_data <- corr_data %>%
     dplyr::select(any_of(c(numeric_cols, ordered_cols)))

if (ncol(corr_data) > 1) {
    corr_data <- corr_data %>%
        mutate(across(all_of(ordered_cols), ~ as.numeric(.)))
    corr_matrix <- stats::cor(corr_data, use = "pairwise.complete.obs")
    corr_pairs <- as.data.frame(as.table(corr_matrix)) %>%
        filter(as.character(Var1) < as.character(Var2)) %>%
        mutate(abs_corr = abs(Freq)) %>%
        arrange(desc(abs_corr))
    high_corr_pairs <- corr_pairs %>%
        filter(abs_corr >= 0.7)
    message("Top correlated predictor pairs:")
    print(utils::head(corr_pairs, 10))
    if (nrow(high_corr_pairs) == 0) {
        message("No predictor pairs exceed |correlation| >= 0.7.")
    } else {
        message("Predictor pairs exceeding |correlation| >= 0.7:")
        print(high_corr_pairs)
    }
} else {
    high_corr_pairs <- data.frame()
    message("Not enough predictor columns to compute correlations.")
}

# Iterate over all response_vars
results_list <- foreach(response_var = response_vars) %do% {
    boot_forward_list <- bootstrap_models(data_df, variables, n_boot, response_var, seed=1337, steps=1000)
    saveRDS(boot_forward_list, file = paste0("discrepancies/2025-10-Data-Version/processed_data/bootstrapped_models_forward_", response_var, ".rds"))
    variable_freq_results <- plot_variable_frequency(boot_forward_list, response_var, data_df, n_boot)
    combination_results <- plot_model_combinations(boot_forward_list, response_var, data_df, n_boot)
    
    list(
        boot_forward_list = boot_forward_list,
        variable_freq = variable_freq_results,
        combination = combination_results
    )

}

stopCluster(cl)

