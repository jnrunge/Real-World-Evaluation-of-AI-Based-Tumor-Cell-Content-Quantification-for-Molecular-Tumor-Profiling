if(file.exists("code/00-universal-dependencies.R")) {
    setwd(file.path(getwd(), "code"))
} 

project_dir <- getwd()

library(renv)

renv::load(".")

# this was decided after reviewing the univariate analyses
# which still contain these variables
# they were then removed because of their low relation with discrepancy 
# and / or their high correlation with other variables
# to improve model selection stability
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

# these were removed from the path vs fmi analyses
# because they a priori do not make sense in that context
variables_not_in_path_vs_fmi <- c("Path_WSI Quality", "AI_scanning_artifacts_area_%")

library(tidyverse)
library(readxl)
library(GGally)
library(patchwork)
library(sjPlot)
library(broom)
library(ggpubr)
library(purrr)
library(foreach)
library(doParallel)
library(vip)
library(ggdendro)
library(dbscan)
library(viridis)
library(factoextra)
library(uwot)
library(ggnewscale)
library(ggforce)
library(parallel)



print(paste0("Project directory is: ", project_dir))

pretty_names <- read_csv(file.path(project_dir, "input/variable_pretty_names.csv"))
# Helper function to clean variable names
clean_var_names <- function(x) {
    x <- stringr::str_replace_all(x, "Foundation.Model...TME.v1.0.1.alpha...Supporting.Result...", "")
    x <- stringr::str_replace_all(x, "Foundation.Model...TME.v1.0.1.alpha...Key.Result...", "")
    x <- stringr::str_replace_all(x, "Artifact Detect v3.0.0 - Key Result - ", "")
    x <- stringr::str_replace_all(x, "TumorDetect v1.2.0 - Supporting Result - ", "")
    x <- stringr::str_replace_all(x, "Artifact Detect v3.0.0 - Supporting Result - ", "")
    return(x)
}

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
    # fallback
    #write_csv(tibble(variable = var, pretty_name = var), "discrepancies/variable_pretty_names.csv", append = TRUE)
    var
}

short_var_label <- function(x) {
    abbrev <- c(
        "inflammation" = "infl.",
        "content" = "cont.",
        "Tertiary Lymphoid Structures" = "TLS",
        "Tumor Infiltrating lymphocytes" = "TILs",
        "Highly cellular stroma" = "cell. stroma"
    )
    x %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_replace_all(" +", " ") %>%
        {
            tmp <- .
            for (pat in names(abbrev)) {
                tmp <- stringr::str_replace_all(
                    tmp,
                    stringr::regex(pat, ignore_case = TRUE), abbrev[pat]
                )
            }
            tmp
        } %>%
        stringr::str_trim()
}

