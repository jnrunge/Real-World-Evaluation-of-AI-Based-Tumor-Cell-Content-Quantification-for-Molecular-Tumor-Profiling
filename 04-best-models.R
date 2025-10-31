library(vip)
library(sjPlot)
pretty_names <- read_csv("discrepancies/variable_pretty_names.csv")
source("discrepancies/2025-10-Data-Version/functions/best_models.R")
response_vars <- variables %>%
    filter(type == "response") %>%
    pull(variable)

classic_no_forced_interactions <- list()

for (response_var in response_vars) {
    var_table <- variables %>%
        filter(type != "response" | variable == response_var)
    
    classic_no_forced_interactions[[response_var]] <- perform_forward_selection(var_table %>% filter(type == "response" | !grepl("TCC", variable)), data_df, force_interaction = FALSE, steps = 5000, seed = 1, reps = 32)
}



## interpret models 
# use compare_model_drop_terms(model, terms) to get term importance/significance