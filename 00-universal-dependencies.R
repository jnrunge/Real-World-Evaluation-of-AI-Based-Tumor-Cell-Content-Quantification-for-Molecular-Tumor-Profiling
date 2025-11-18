################################################################################
# Script: 00-universal-dependencies.R
# Purpose: Load universal dependencies, helper functions, and variable mappings
#          for TCC discrepancy analysis
# Author: Jan-Niklas Runge
# 
# Description:
#   - Sets up project directory structure
#   - Loads required R packages via renv
#   - Defines excluded variables for modeling
#   - Provides helper functions for variable name formatting
# 
# Inputs:
#   - input/variable_pretty_names.csv: Mapping of technical to display names
# 
# Outputs:
#   - Loaded libraries and helper functions in global environment
#   - project_dir: Path to project root
# 
# Dependencies: See renv.lock for package versions
################################################################################

# Set working directory - adjust if running from different location
if(file.exists("code/00-universal-dependencies.R")) {
    setwd(file.path(getwd(), "code"))
} 

project_dir <- getwd()

# Load renv environment for reproducibility
library(renv)
renv::load(".")

# Variables excluded from modeling ---------------------------------------------
# Exclusion rationale: Low association with discrepancy outcomes and/or
# high correlation with other predictors (reduces model selection stability)
# Decision made after reviewing univariate analyses
variables_not_in_model <- c("Metastatic",
"SpecimenSiteGrouped_01",
"SpecimenSiteGrouped_02",
"AI_total_tumor_area_%",
"Path_Cautery/Crush artifact average_yes_no",
"Path_Tumor-Stroma Content",
"Path_High content of intratumoral immune infiltrates",
"Path_Highly cellular stroma",
"AI_cancer_area_%_per_tumor",
"AI_stroma_%_per_tumor",
"AI_lymphocyte_density_in_tumor_mm2",
"AI_fibroblast_density_in_tumor_mm2",
"AI_%_lymphocytes_in_tumor",
"AI_%_of_fibroblasts_in_tumor",
"AI_%_of_immune_cells_in_tumor",
"AI_%_of_non_immune_cells_in_tumor",
"Path_Margin Ink_Resection"
)

# Variables excluded from Path vs FMI analyses ---------------------------------
# Exclusion rationale: Image quality metrics not applicable to FMI comparisons
variables_not_in_path_vs_fmi <- c("Path_WSI Quality", "AI_scanning_artifacts_area_%")

# Load required packages -------------------------------------------------------
required_packages <- c(
    "tidyverse", "readxl", "GGally", "patchwork", "sjPlot", "broom",
    "ggpubr", "purrr", "foreach", "doParallel", "vip", "ggdendro",
    "dbscan", "viridis", "factoextra", "uwot", "ggnewscale", 
    "ggforce", "parallel"
)

for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        stop("Required package '", pkg, "' is not installed. ",
             "Please run renv::restore() to install all dependencies.")
    }
}

print(paste0("Project directory is: ", project_dir))

# Load variable name mappings --------------------------------------------------
pretty_names_path <- file.path(project_dir, "input/variable_pretty_names.csv")
if (!file.exists(pretty_names_path)) {
    stop("Required file not found: ", pretty_names_path, "\n",
         "Please ensure input/variable_pretty_names.csv exists in project directory.")
}
pretty_names <- read_csv(pretty_names_path, show_col_types = FALSE)

# Helper function to clean variable names
clean_var_names <- function(x) {
    x <- stringr::str_replace_all(x, "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...", "")
    x <- stringr::str_replace_all(x, "Foundation.Model...TME.v1.0.1.alpha...Key.Result...", "")
    x <- stringr::str_replace_all(x, "Artifact Detect v3.0.0 - Key Result - ", "")
    x <- stringr::str_replace_all(x, "TumorDetect v1.2.0 - Supporting Result - ", "")
    x <- stringr::str_replace_all(x, "Artifact Detect v3.0.0 - Supporting Result - ", "")
    return(x)
}

#' Convert variable names to pretty display names
#' 
#' @param var Character string, variable name (supports interaction terms with ":")
#' @return Character string, formatted display name
#' @details Handles interaction terms by splitting on ":" and formatting each part.
#'          Uses exact matching or prefix matching against pretty_names lookup table.
get_pretty_name <- function(var) {
    # Handle interaction terms like "Var1:Var2"
    if (grepl(":", var, fixed = TRUE)) {
        parts <- strsplit(var, ":", fixed = TRUE)[[1]]
        pretty_parts <- vapply(parts, get_pretty_name, character(1))
        return(paste(pretty_parts, collapse = ":"))
    }
    # Exact match
    if (var %in% pretty_names$variable) {
        return(pretty_names$pretty_name[match(var, pretty_names$variable)])
    }
    # Try to match a prefix in pretty_names$variable
    candidates <- pretty_names$variable[startsWith(var, pretty_names$variable)]
    if (length(candidates) > 0) {
        # pick the longest matching prefix
        best <- candidates[which.max(nchar(candidates))]
        pretty_prefix <- pretty_names$pretty_name[pretty_names$variable == best]
        suffix <- substring(var, nchar(best) + 1)
        # strip leading separators from suffix
        suffix <- sub("^[ _:\\.]+", "", suffix)
        return(paste0(pretty_prefix, if (nzchar(suffix)) suffix else ""))
    }
    # fallback: return original variable name
    var
}


# Print session info for reproducibility ---------------------------------------
message("\n=== Session Info ===")
print(sessionInfo())

