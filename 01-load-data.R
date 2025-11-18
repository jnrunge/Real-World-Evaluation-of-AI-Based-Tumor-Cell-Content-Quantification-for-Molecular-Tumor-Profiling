################################################################################
# Script: 01-load-data.R
# Purpose: Load and preprocess PathAI dataset for TCC discrepancy analysis
# Author: Jan-Niklas Runge
# 
# Description:
#   - Loads raw PathAI dataset (300 cases) and legacy merged data
#   - Processes variables: converts types, orders factors, handles interactions
#   - Identifies outliers but does not remove them
#   - Scales numeric variables for modeling
#   - Generates diagnostic plots (distributions, outliers)
# 
# Inputs:
#   - input/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx (sheets 1-2)
#   - input/raw/merged_data.csv (rawer data for AI artifact measures)
# 
# Outputs:
#   - output/processed_data/data_df.rds: Scaled, cleaned dataset
#   - output/processed_data/data_df_pre_scaling.rds: Pre-scaling version
#   - output/processed_data/data_df_renamed.rds: R-safe column names
#   - output/processed_data/data_df_complete.rds: Full dataset
#   - output/processed_data/data_df_pre_scaling_NAd_sampletypes.rds: With NAs for variables out of sample type scope
#   - output/processed_data/variables.rds: Variable lists for modeling
#   - output/outliers/*.png: Outlier diagnostic plots
#   - output/distributions/*.png: Distribution plots for all variables
# 
# Dependencies: See 00-universal-dependencies.R for package requirements
################################################################################

# Configuration ----------------------------------------------------------------
# IQR-based outlier detection bounds (percentiles and multiplier)
OUTLIER_LOWER_PERCENTILE <- 0.05
OUTLIER_UPPER_PERCENTILE <- 0.95
OUTLIER_IQR_MULTIPLIER <- 1

# Minimum unique values for outlier detection (excludes sparse variables)
MIN_UNIQUE_FOR_OUTLIER_DETECTION <- 10

# Check if cached output files exist -------------------------------------------
output_files_exist <- all(file.exists(
  file.path(project_dir, "output/processed_data/data_df_pre_scaling.rds"),
  file.path(project_dir, "output/processed_data/data_df.rds"),
  file.path(project_dir, "output/processed_data/data_df_renamed.rds"),
  file.path(project_dir, "output/processed_data/variables.rds"),
  file.path(project_dir, "output/processed_data/data_df_complete.rds"),
  file.path(project_dir, "output/processed_data/data_df_pre_scaling_NAd_sampletypes.rds")
))

# Load variable metadata -------------------------------------------------------
var_desc <- readxl::read_excel(
  file.path(project_dir, "input/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx"), 
  sheet = 2
)

# Helper functions -------------------------------------------------------------

#' Determine relevant sample type categories for a variable
#' 
#' @param var_name Character string, variable name
#' @return Character vector of relevant sample types, or FALSE if not applicable
#' @details Special handling for cancer content ratio; otherwise uses metadata
get_relevant_categories <- function(var_name) {
  # Special case: cancer content ratio applies to both Biopsy and Resection
  if (var_name %in% c("Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed")) {
    return(c("Biopsy", "Resection"))
  }
  
  if (get_interaction_or_not(var_name)) {
    sample_types <- var_desc %>% 
      filter(`Independent Preanalytical Variables` == var_name) %>% 
      pull(`Sample type relevance`)
    return(strsplit(sample_types, "/")[[1]] %>% str_trim())
  } else {
    return(FALSE)
  }
}

#' Check if variable has sample-type specific interactions
#' 
#' @param var_name Character string, variable name
#' @return Logical, TRUE if variable is sample-type specific
#' @details Returns FALSE for variables marked "All" or NA in metadata
get_interaction_or_not <- function(var_name) {
  # Special case variables
  if (var_name %in% c("Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed")) {
    return(TRUE)
  }
  
  # Check metadata
  if (var_name %in% (var_desc[, 2] %>% pull())) {
    relevance <- var_desc %>% 
      filter(`Independent Preanalytical Variables` == var_name) %>% 
      pull(`Sample type relevance`)
    return(!any(c(NA, "", "All") %in% relevance))
  } else {
    return(FALSE)
  }
}

# Load or process data ---------------------------------------------------------

if (output_files_exist) {
  message("Output files already exist. Loading from saved RDS files...")
  data_df_pre_scaling <- readRDS(file.path(project_dir, "output/processed_data/data_df_pre_scaling.rds"))
  data_df <- readRDS(file.path(project_dir, "output/processed_data/data_df.rds"))
  data_df_renamed <- readRDS(file.path(project_dir, "output/processed_data/data_df_renamed.rds"))
  data_df_complete <- readRDS(file.path(project_dir, "output/processed_data/data_df_complete.rds"))
  variables <- readRDS(file.path(project_dir, "output/processed_data/variables.rds"))
  data_df_pre_scaling_NAd_sampletypes <- readRDS(file.path(project_dir, "output/processed_data/data_df_pre_scaling_NAd_sampletypes.rds"))
  message("Data loaded successfully.")
} else {
  message("Output files not found. Running full data processing pipeline...")

  # Load raw data --------------------------------------------------------------
  data_df <- read_excel(
    file.path(project_dir, "input/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx"), 
    sheet = 1
  )
  
  # Load legacy data for artifact measures
  old_data <- read_csv2(file.path(project_dir, "input/raw/merged_data.csv"))

  # Compute derived variables --------------------------------------------------
  
  # Metastatic status: cancer site differs from specimen site
  data_df$Metastatic <- data_df$Cancer != data_df$SpecimenSite
  
  # Fix naming: use non-rounded AI tumor cell content
  data_df$TCC_AI <- data_df$TC_AI
  
  # Compute discrepancy variables (primary outcomes)
  data_df$TCC_Patho_minus_TCC_AI <- data_df$TCC_Patho - data_df$TCC_AI
  data_df$TCC_FMI_minus_TCC_AI <- data_df$TCC_FMI - data_df$TCC_AI

  # Add artifact fraction measure from AI artifact detection
  old_data$AI_artefact_fraction_measure <- (
    1 - (old_data$`Artifact Detect v3.0.0 - Supporting Result - Area of Evaluable Tissue - mm^2` / 
         old_data$`Artifact Detect v3.0.0 - Supporting Result - Area of Total Tissue - mm^2`)
  )
  data_df <- left_join(
    data_df, 
    old_data %>% select(pseudonym, AI_artefact_fraction_measure), 
    by = c("Pseudonym" = "pseudonym")
  )

  # Fix fragmented cancer content ratio variable
  # Rationale: Original variable has "fragmented" entries; recalculate from components
  old_data$`Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed` <- 
    old_data$`Cancer content N specimens_Resection/Biopsy_numeric` / 
    old_data$`N Different tissue pieces per slide_Resection/Biopsy_numeric`
  old_data$`Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed`[
    is.na(old_data$`Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed`)
  ] <- 0
  
  data_df <- left_join(
    data_df, 
    old_data %>% select(pseudonym, `Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed`), 
    by = c("Pseudonym" = "pseudonym")
  )

  # Verify data integrity against legacy dataset ------------------------------
  # Compare columns between data_df and old_data to identify matching ones
  joined_data <- full_join(data_df, old_data, by = c("Pseudonym" = "pseudonym"))

  common_cols <- map_dfr(names(data_df), function(col1) {
    if (col1 == "Pseudonym") return(NULL)

    map_dfr(names(old_data), function(col2) {
      if (col2 == "pseudonym") return(NULL)

      # Get values from joined data (already aligned by Pseudonym)
      vec1 <- joined_data[[if (paste0(col1, ".x") %in% names(joined_data)) paste0(col1, ".x") else col1]]
      vec2 <- joined_data[[if (paste0(col2, ".y") %in% names(joined_data)) paste0(col2, ".y") else col2]]

      # Remove NAs from both
      valid_idx <- !is.na(vec1) & !is.na(vec2)
      vec1 <- vec1[valid_idx]
      vec2 <- vec2[valid_idx]

      if (length(vec1) > 0 && identical(vec1, vec2)) {
        tibble(data_df_col = col1, old_data_col = col2)
      }
    })
  })

  if (nrow(common_cols) > 0) {
    message("Matching columns found between datasets:")
    print(common_cols)
  } else {
    message("No matching columns found between datasets")
  }

  # Extract independent variables from metadata -------------------------------
  list_of_all_indepedent_variables <- read_excel(
    file.path(project_dir, "input/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx"), 
    sheet = 2
  ) %>% select(2)
  
  list_of_all_indepedent_variables <- pull(list_of_all_indepedent_variables)
  
  # Extract variables up to "Dependent variables" marker
  list_of_all_indepedent_variables <- list_of_all_indepedent_variables[
    1:(which(list_of_all_indepedent_variables == "Dependent variables") - 1)
  ]
  
  # Remove NA and section headers
  list_of_all_indepedent_variables <- list_of_all_indepedent_variables[!is.na(list_of_all_indepedent_variables)]
  list_of_all_indepedent_variables <- list_of_all_indepedent_variables[
    !list_of_all_indepedent_variables %in% c("Independent Pre-/Analytical variables", "Independent Analytical Variables")
  ]
  
  # Add custom variables
  list_of_all_indepedent_variables <- c(list_of_all_indepedent_variables, "AI_artefact_fraction_measure")
  
  # Replace old variable name with fixed version
  list_of_all_indepedent_variables <- str_replace(
    list_of_all_indepedent_variables,
    "Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio",
    "Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed"
  )

  # Validate all required variables exist in dataset
  missing_vars <- list_of_all_indepedent_variables[
    !list_of_all_indepedent_variables %in% colnames(data_df)
  ]
  
  if (length(missing_vars) > 0) {
    stop("The following variables are missing from the data: ", paste0(missing_vars, collapse = ", "))
  }

  # Create variable metadata table ---------------------------------------------
  variables <- tibble(variable = list_of_all_indepedent_variables, type = "explanatory") %>%
    rowwise() %>%
    mutate(interaction = if_else(get_interaction_or_not(variable), "Sample type", "")) %>%
    ungroup()
  
  # Add response variables (discrepancy measures)
  variables <- bind_rows(
    variables, 
    tibble(
      variable = c("TCC_Patho_minus_TCC_AI", "TCC_FMI_minus_TCC_AI", "TCC_Patho_minus_TCC_FMI"), 
      type = "response"
    )
  )

  # Rename and subset data -----------------------------------------------------
  data_df <- data_df %>% rename(TCC_Patho_minus_TCC_FMI = TC_Patho_minus_TCC_FMI)
  data_df_complete <- data_df
  data_df <- data_df %>% select(all_of(variables$variable))

  # Convert numeric-like character columns to numeric -------------------------
  # Rationale: Some columns stored as character can be safely converted to numeric
  numeric_like_columns <- names(data_df)[sapply(data_df, function(v) {
    if (is.numeric(v)) return(FALSE)
    
    # Try for comma-separated decimal strings (European format)
    if (all(str_detect(v, "[0-9]+(,)[0-9]+"))) {
      v <- str_replace_all(v, ",", ".")
    }
    
    suppressWarnings(v_num <- as.numeric(v))
    na_orig <- sum(is.na(v))
    na_num <- sum(is.na(v_num))
    n_unique <- length(unique(v_num[!is.na(v_num)]))
    total <- length(v)
    
    # Convert if <1% conversion loss and >5 unique values
    (na_num - na_orig) / total < 0.01 && n_unique > 5
  })]

  if (length(numeric_like_columns) > 0) {
    message("Converting character columns to numeric:")
    print(numeric_like_columns)
    
    data_df <- data_df %>%
      mutate(across(all_of(numeric_like_columns), ~ {
        # Handle European decimal format (comma as separator)
        if (all(str_detect(na.omit(.), "[0-9]+(,)[0-9]+"))) {
          as.numeric(str_replace_all(., ",", "."))
        } else {
          as.numeric(.)
        }
      }))
  }

  # Order categorical variables ------------------------------------------------
  # Define ordinal level orderings for common factor types
  quality_order <- c("Poor", "Satisfactory", "Good")
  absent_present_order <- c("Absent", "Minimal", "Prominent")
  no_yes_order <- c("No", "Yes")
  lowhigh_order <- c("Low", "High")
  lowhigh_order_2 <- c("No", "Low", "Intermediate", "High")
  absent_present_order_2 <- c("Absent", "Minimal", "Moderate", "Extensive")

  orders_list <- list(
    quality_order, 
    absent_present_order, 
    absent_present_order_2, 
    no_yes_order, 
    lowhigh_order, 
    lowhigh_order_2
  )

  # Apply orderings to character columns
  for (var in names(data_df)[sapply(data_df, is.character)]) {
    vals <- na.omit(unique(data_df[[var]]))
    
    # Try each ordering; once matched, convert and move on
    for (ord in orders_list) {
      if (length(vals) > 0 && all(vals %in% ord)) {
        data_df[[var]] <- factor(data_df[[var]], levels = ord, ordered = TRUE)
        break
      }
    }
  }

  # Check for unordered categorical variables
  unordered_vars <- names(data_df)[sapply(data_df, function(v) {
    (is.character(v) || is.factor(v)) && !is.ordered(v)
  })]

  if (length(unordered_vars) > 0) {
    message("Unordered character/factor variables:")
    print(unordered_vars)
  } else {
    message("No unordered character/factor variables found")
  }

  # Handle sample-type specific variables --------------------------------------
  # Rationale: Some variables only apply to specific sample types (Biopsy/Resection)
  # Set out-of-scope values to baseline (first level for factors, 0 for numeric)
  
  relevant_categories <- lapply(colnames(data_df), get_relevant_categories)
  names(relevant_categories) <- colnames(data_df)
  relevant_categories <- Filter(function(x) !identical(x, FALSE), relevant_categories)
  
  types_tbl <- tibble(variable = names(relevant_categories)) %>%
    mutate(
      typeof = map_chr(variable, ~ typeof(data_df[[.x]])),
      class  = map_chr(variable, ~ paste(class(data_df[[.x]]), collapse = "/"))
    )

  print(types_tbl, n = nrow(types_tbl))

  # For all variables with sample-type relevance, set out-of-scope rows to baseline
  for (i in seq_len(nrow(types_tbl))) {
    var <- types_tbl$variable[i]
    cls <- types_tbl$class[i]
    rel <- relevant_categories[[var]]
    if (is.null(rel)) next

    idx <- !(data_df$`Sample type` %in% rel)

    if (cls == "ordered/factor") {
      # Assign lowest level
      lvl1 <- levels(data_df[[var]])[1]
      data_df[[var]][idx] <- lvl1
    } else if (cls == "numeric") {
      data_df[[var]][idx] <- 0
    }
  }

  # Create version with sample-type NAs (for diagnostics) ---------------------
  # Rationale: Preserve a version where out-of-scope values are NA, not baseline
  data_df_pre_scaling_NAd_sampletypes <- data_df %>%
    mutate(across(everything(), ~ {
      relevant_categories <- get_relevant_categories(cur_column())
      if (isFALSE(relevant_categories)) return(.)
      
      # Preserve the original class/type
      original_col <- .
      mask <- `Sample type` %in% relevant_categories
      original_col[!mask] <- NA
      original_col
    }))

  # Remove rows with NA values -------------------------------------------------
  # Rationale: Complete case analysis required for modeling
  rows_with_na <- which(apply(data_df, 1, function(row) any(is.na(row))))

  if (length(rows_with_na) > 0) {
    message(paste("Found", length(rows_with_na), "rows with NA values"))
    message("Row indices with NAs: ", paste(rows_with_na, collapse = ", "))

    # Show which variables have NAs in each row
    for (row_idx in rows_with_na) {
      na_vars <- names(data_df)[is.na(data_df[row_idx, ])]
      message(paste("  Row", row_idx, "has NAs in:", paste(na_vars, collapse = ", ")))
    }

    # Remove rows with NAs
    data_df <- data_df[!apply(data_df, 1, function(row) any(is.na(row))), ]
    message(paste("Removed", length(rows_with_na), "rows. Remaining rows:", nrow(data_df)))
  } else {
    message("No rows with NA values found")
  }

  # Identify outliers using IQR method -----------------------------------------
  numeric_vars <- names(data_df)[sapply(data_df, is.numeric)]

  # For each numeric variable, compute IQR bounds and collect outlier row indices
  outliers <- map_dfr(numeric_vars, function(var) {
    x <- data_df_pre_scaling_NAd_sampletypes[[var]]
    
    # Skip too-sparse or constant vectors
    if (length(unique(na.omit(x))) < MIN_UNIQUE_FOR_OUTLIER_DETECTION) {
      return(NULL)
    }
    
    qnt <- quantile(x, c(OUTLIER_LOWER_PERCENTILE, OUTLIER_UPPER_PERCENTILE), na.rm = TRUE)
    iqr <- qnt[2] - qnt[1]
    lb <- qnt[1] - OUTLIER_IQR_MULTIPLIER * iqr
    ub <- qnt[2] + OUTLIER_IQR_MULTIPLIER * iqr
    idx <- which(!is.na(x) & (x < lb | x > ub))
    
    if (length(idx) == 0) return(NULL)
    
    tibble(
      variable    = var,
      row         = idx,
      value       = x[idx],
      lower_bound = lb,
      upper_bound = ub
    )
  })

  # Summary: how many outliers per variable
  message("\nOutlier summary:")
  outliers %>%
    count(variable, name = "n_outliers") %>%
    mutate(pct = n_outliers / nrow(data_df_pre_scaling_NAd_sampletypes) * 100) %>%
    arrange(desc(n_outliers)) %>%
    print()

  outlier_rows <- unique(outliers$row)

  # Generate outlier diagnostic plots ------------------------------------------
  output_dir <- file.path(project_dir, "output/outliers/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (var in unique(outliers$variable)) {
    df <- tibble(
      row = seq_len(nrow(data_df_pre_scaling_NAd_sampletypes)),
      value = data_df_pre_scaling_NAd_sampletypes[[var]]
    )
    df_out <- outliers %>% filter(variable == var)
    
    p <- ggplot(df, aes(x = row, y = value)) +
      geom_point(alpha = 0.6) +
      geom_point(
        data = df_out, aes(x = row, y = value),
        colour = "red", size = 2
      ) +
      labs(
        title = paste("Outliers in", var),
        x = "Row index", 
        y = "Value"
      ) +
      theme_bw(14)
    
    ggsave(
      filename = file.path(output_dir, paste0(gsub("[^A-Za-z0-9]", "_", var), ".png")), 
      plot = p, 
      width = 12, 
      height = 4
    )
  }

  # Generate distribution plots ------------------------------------------------
  output_dir <- file.path(project_dir, "output/distributions/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (var in names(data_df_pre_scaling_NAd_sampletypes)) {
    p <- ggplot(data_df_pre_scaling_NAd_sampletypes, aes_string(x = paste0("`", var, "`"))) +
      theme_bw(14) +
      labs(title = paste("Distribution of", var))

    if (is.numeric(data_df_pre_scaling_NAd_sampletypes[[var]])) {
      p <- p + geom_histogram(bins = 50, fill = "blue", color = "black", alpha = 0.7)
    } else if (is.factor(data_df_pre_scaling_NAd_sampletypes[[var]]) || 
               is.character(data_df_pre_scaling_NAd_sampletypes[[var]])) {
      p <- p + geom_bar(fill = "blue", color = "black", alpha = 0.7)
    } else {
      next # Skip if neither numeric nor categorical
    }

    ggsave(
      filename = file.path(output_dir, paste0(gsub("[^A-Za-z0-9]", "_", var), ".png")), 
      plot = p, 
      width = 13, 
      height = 6
    )
  }
  
  # Final validation -----------------------------------------------------------
  na_columns <- colnames(data_df)[sapply(data_df, function(x) any(is.na(x)))]
  if (length(na_columns) > 0) {
    message("Columns with NAs found:")
    print(na_columns)
  } else {
    message("No columns with NAs found.")
  }

  # Scale numeric variables ----------------------------------------------------
  data_df_pre_scaling <- data_df
  
  # Save pre-scaling versions
  dir.create(file.path(project_dir, "output/processed_data/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(data_df_pre_scaling_NAd_sampletypes, 
          file = file.path(project_dir, "output/processed_data/data_df_pre_scaling_NAd_sampletypes.rds"))
  saveRDS(data_df_pre_scaling, 
          file = file.path(project_dir, "output/processed_data/data_df_pre_scaling.rds"))
  
  # Scaling function: center and scale to unit variance
  scale_ <- function(x) {
    return(scale(x, center = TRUE, scale = TRUE) %>% as.numeric())
  }

  response_vars <- variables %>%
    filter(type == "response") %>%
    pull(variable)

  # Scale all numeric variables except response variables
  data_df <- data_df %>%
    mutate(across(all_of(setdiff(numeric_vars, response_vars)), scale_))

  # Optional: Remove highly correlated/redundant variables --------------------
  # Uncomment to exclude highly correlated predictors
  # highly_cor_vars <- c("AI_non_immune_cell_density_in_tumor_mm2",
  #                      "AI_immune_cell_density_in_tumor_mm2", 
  #                      "AI_%_lymphocytes_in_tumor",
  #                      "AI_%_of_fibroblasts_in_tumor", 
  #                      "AI_%_of_immune_cells_in_tumor", 
  #                      "AI_%_of_non_immune_cells_in_tumor")
  # variables <- variables %>% filter(!variable %in% highly_cor_vars)
  # data_df <- data_df %>% dplyr::select(-all_of(highly_cor_vars))

  # Save processed datasets ----------------------------------------------------
  data_df_renamed <- data_df
  names(data_df_renamed) <- make.names(names(data_df_renamed), unique = TRUE)
  
  saveRDS(data_df_complete, 
          file = file.path(project_dir, "output/processed_data/data_df_complete.rds"))
  saveRDS(data_df, 
          file = file.path(project_dir, "output/processed_data/data_df.rds"))
  saveRDS(data_df_renamed, 
          file = file.path(project_dir, "output/processed_data/data_df_renamed.rds"))
  saveRDS(variables, 
          file = file.path(project_dir, "output/processed_data/variables.rds"))

  message("Data processing complete and saved.")
}

# Extract response variable names for downstream use --------------------------
response_vars <- variables %>%
  filter(type == "response") %>%
  pull(variable)
