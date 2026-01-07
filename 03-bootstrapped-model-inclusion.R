################################################################################
# Script: 03-bootstrapped-model-inclusion.R
# Purpose: Perform bootstrapped stepwise regression for variable selection stability
# Author: Jan-Niklas Runge
# 
# Description:
#   - Runs forward stepwise selection on n_boot bootstrap samples
#   - Quantifies variable inclusion frequency across bootstrap iterations
# 
# 
# Outputs:
#   - output/processed_data/bootstrapped_models_forward_[outcome].rds: Bootstrap results
#   - output/bootstrapped_models/variable_frequency_[outcome].pdf: Inclusion frequency plots
# 
################################################################################

# Configuration ----------------------------------------------------------------
# Number of stepwise iterations per bootstrap sample
# Rationale: Ensures convergence of forward selection algorithm
STEPS_TO_RUN <- 5000

# Number of bootstrap samples
# Rationale: 1000 iterations provides stable frequency estimates
N_BOOT <- 1000

# Correlation threshold for multicollinearity detection
# Rationale: |r| >= 0.7 indicates potential redundancy between predictors
HIGH_CORR_THRESHOLD <- 0.7

# Random seed for reproducibility
RANDOM_SEED <- 1337

# Parallel processing: reserve 1 core for system
N_CORES <- parallel::detectCores() - 1

# Load helper functions --------------------------------------------------------
source(file.path(project_dir, "functions/manual-stepwise.R"))
source(file.path(project_dir, "functions/bootstrapping_models.R"))

# Extract response variables ---------------------------------------------------
response_vars <- variables %>%
    filter(type == "response") %>%
    pull(variable)

# Apply exclusions to predictor set -------------------------------------------
# Rationale: Remove variables excluded based on univariate analyses
variables <- variables %>%
    filter(!(variable %in% variables_not_in_model))

# Set up parallel processing ---------------------------------------------------
message(paste("Initializing parallel processing with", N_CORES, "cores"))
cl <- makeCluster(N_CORES)
registerDoParallel(cl)

# Detect highly correlated predictors -----------------------------------------
# Rationale: High correlations can destabilize stepwise selection and inflate VIF
message("\n=== Checking for highly correlated predictor pairs ===")

predictor_vars <- variables %>%
    filter(type != "response") %>%
    pull(variable)

available_predictors <- intersect(predictor_vars, names(data_df))

# Select numeric and ordered factors for correlation analysis
corr_data <- data_df %>%
    dplyr::select(any_of(available_predictors))

ordered_cols <- names(corr_data)[vapply(corr_data, is.ordered, logical(1))]
numeric_cols <- names(corr_data)[vapply(corr_data, is.numeric, logical(1))]
corr_data <- corr_data %>%
    dplyr::select(any_of(c(numeric_cols, ordered_cols)))

if (ncol(corr_data) > 1) {
    # Convert ordered factors to numeric for correlation computation
    corr_data <- corr_data %>%
        mutate(across(all_of(ordered_cols), ~ as.numeric(.)))
    
    # Compute pairwise correlations
    corr_matrix <- stats::cor(corr_data, use = "pairwise.complete.obs")
    
    # Extract upper triangle (avoid duplicates)
    corr_pairs <- as.data.frame(as.table(corr_matrix)) %>%
        filter(as.character(Var1) < as.character(Var2)) %>%
        mutate(abs_corr = abs(Freq)) %>%
        arrange(desc(abs_corr))
    
    # Identify high-correlation pairs
    high_corr_pairs <- corr_pairs %>%
        filter(abs_corr >= HIGH_CORR_THRESHOLD)
    
    message("Top correlated predictor pairs:")
    print(utils::head(corr_pairs, 10))
    
    if (nrow(high_corr_pairs) == 0) {
        message(paste("No predictor pairs exceed |correlation| >=", HIGH_CORR_THRESHOLD))
    } else {
        message(paste("Predictor pairs exceeding |correlation| >=", HIGH_CORR_THRESHOLD, ":"))
        print(high_corr_pairs)
    }
} else {
    high_corr_pairs <- data.frame()
    message("Not enough predictor columns to compute correlations.")
}

# Create output directory ------------------------------------------------------
output_dir <- file.path(project_dir, "output/bootstrapped_models")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Run bootstrap analysis for each outcome -------------------------------------
message("\n=== Starting bootstrap model selection ===")

results_list <- foreach(response_var = response_vars) %do% {
    message(paste("\n--- Processing outcome:", response_var, "---"))
    
    # Apply Path vs FMI exclusions if applicable
    # Rationale: Image quality variables irrelevant for FMI comparisons
    current_variables <- variables
    if (response_var == "TCC_Patho_minus_TCC_FMI") {
        current_variables <- current_variables %>%
            filter(!(variable %in% variables_not_in_path_vs_fmi))
        message("Applied Path vs FMI exclusions")
    }
    
    # Define output file path
    rds_file <- file.path(project_dir, "output/processed_data", 
                         paste0("bootstrapped_models_forward_", response_var, ".rds"))
    
    # Check for existing results -----------------------------------------------
    if (file.exists(rds_file)) {
        boot_forward_list <- readRDS(rds_file)
        
        # Validate number of bootstrap iterations
        if (length(boot_forward_list) != N_BOOT) {
            message(paste("Number of bootstraps mismatch for", response_var, 
                         ". Expected:", N_BOOT, "Found:", length(boot_forward_list)))
            message("Rerunning bootstrap.")
            boot_forward_list <- bootstrap_models(
                data_df, current_variables, N_BOOT, response_var, 
                seed = RANDOM_SEED, steps = STEPS_TO_RUN
            )
            saveRDS(boot_forward_list, file = rds_file)
        } else {
            message(paste("Loading existing bootstrap results for", response_var))
        }
    } else {
        message(paste("No existing bootstrap results for", response_var, ". Running bootstrap."))
        boot_forward_list <- bootstrap_models(
            data_df, current_variables, N_BOOT, response_var, 
            seed = RANDOM_SEED, steps = STEPS_TO_RUN
        )
        saveRDS(boot_forward_list, file = rds_file)
    }
    
    # Generate variable frequency plots ----------------------------------------
    # Rationale: Visualize how often each variable is selected across bootstraps
    variable_freq_results <- plot_variable_frequency(
        boot_forward_list, response_var, data_df, N_BOOT
    )
    
    # Model combination analysis (currently disabled) --------------------------
    # Rationale: Analyzing specific predictor combinations can identify 
    # co-selected variable sets, but not currently of interest
    # combination_results <- plot_model_combinations(boot_forward_list, response_var, data_df, N_BOOT)
    
    message(paste("Completed bootstrap analysis for", response_var))
    
    # Return results
    list(
        boot_forward_list = boot_forward_list,
        variable_freq = variable_freq_results
        # combination = combination_results  # Commented out as not used
    )
}

# Clean up parallel processing -------------------------------------------------
stopCluster(cl)

message("\n=== Bootstrap model selection complete ===")
message(paste("Results saved to:", output_dir))

