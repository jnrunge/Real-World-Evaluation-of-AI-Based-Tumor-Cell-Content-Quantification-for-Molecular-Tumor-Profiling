library(tidyverse)
library(readxl)

## The actual df
data_df <- read_excel("data/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx", sheet = 1)
data_df$Metastatic <- data_df$Cancer != data_df$SpecimenSite
old_data <- read_csv2("data/merged_data.csv")

# un-round the TCC_AI variables
data_PathAI <- read_csv("data/20240710_500_samples_with_results.csv")
data_df <- left_join(data_df, data_PathAI %>% select(Pseudonym = pseudonym...1, TCC_AI = Foundation.Model...TME.v1.0.1.alpha...Key.Result...Percent.tumor.nuclei....))
data_df$TCC_Patho_minus_TCC_AI <- data_df$TCC_Patho - data_df$TCC_AI
data_df$TCC_FMI_minus_TCC_AI <- data_df$TCC_FMI - data_df$TCC_AI

# add a general artefact % measure
old_data$AI_artefact_fraction_measure <- (1 - (old_data$`Artifact Detect v3.0.0 - Supporting Result - Area of Evaluable Tissue - mm^2` / old_data$`Artifact Detect v3.0.0 - Supporting Result - Area of Total Tissue - mm^2`))
data_df <- left_join(data_df, old_data %>% select(pseudonym, AI_artefact_fraction_measure), by = c("Pseudonym" = "pseudonym"))

# fix "Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio" variable (has "fragmented" entries)
old_data$`Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed` <- old_data$`Cancer content N specimens_Resection/Biopsy_numeric` / old_data$`N Different tissue pieces per slide_Resection/Biopsy_numeric`
old_data$`Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed`[is.na(old_data$`Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed`)] <- 0
summary(old_data$`Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed`)
data_df <- left_join(data_df, old_data %>% select(pseudonym, `Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed`), by = c("Pseudonym" = "pseudonym"))


## test if anything is different from the raw data (this df was manually handled)
# Compare columns between data_df and old_data to find matching ones
# Join datasets by Pseudonym to align rows
joined_data <- full_join(
  data_df,
  old_data,
  by = c("Pseudonym" = "pseudonym")
)

# Compare columns to find matching ones
common_cols <- map_dfr(names(data_df), function(col1) {
  if (col1 == "Pseudonym") return(NULL)
  
  map_dfr(names(old_data), function(col2) {
    if (col2 == "pseudonym") return(NULL)
    
    # Get values from joined data (already aligned by Pseudonym)
    vec1 <- joined_data[[if(paste0(col1, ".x") %in% names(joined_data)) paste0(col1, ".x") else col1]]
    vec2 <- joined_data[[if(paste0(col2, ".y") %in% names(joined_data)) paste0(col2, ".y") else col2]]
    
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
  message("Matching columns found:")
  print(common_cols)
} else {
  message("No matching columns found")
}

#names(data_df)[!(names(data_df) %in% common_cols$data_df_col)]

## get the list of all desired independent variables
list_of_all_indepedent_variables <- read_excel("data/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx", sheet = 2) %>% select(2)
list_of_all_indepedent_variables <- pull(list_of_all_indepedent_variables)
list_of_all_indepedent_variables <- list_of_all_indepedent_variables[1:(which(list_of_all_indepedent_variables == "Dependent variables") - 1)]
list_of_all_indepedent_variables <- list_of_all_indepedent_variables[!is.na(list_of_all_indepedent_variables)]
list_of_all_indepedent_variables <- list_of_all_indepedent_variables[!list_of_all_indepedent_variables %in% c("Independent Pre-/Analytical variables", "Independent Analytical Variables")]
list_of_all_indepedent_variables <- c(list_of_all_indepedent_variables, "TCC_Patho", "TCC_AI", "TCC_FMI", "AI_artefact_fraction_measure")
# Replace the old variable name with the fixed version
list_of_all_indepedent_variables <- str_replace(
  list_of_all_indepedent_variables,
  "Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio",
  "Path_Cancer_content_in_tissue_specimens_per_slide_Resection/Biopsy_Ratio_fixed"
)

if((missing_vars <- list_of_all_indepedent_variables[which(!list_of_all_indepedent_variables %in% (data_df %>% colnames()))])  %>% length() > 0) {
  stop(paste0("The following variables are missing from the data: ", paste0(missing_vars, collapse = ", ")))
}

variables <- tibble(variable = list_of_all_indepedent_variables, type = "explanatory")
variables <- bind_rows(variables, tibble(variable = c("TCC_Patho_minus_TCC_AI", "TCC_FMI_minus_TCC_AI", "TC_Patho_minus_TCC_FMI"), type = "response"))

data_df <- data_df %>% select(all_of(variables$variable))

numeric_like_columns <- names(data_df)[sapply(data_df, function(v) {
  if (is.numeric(v)) {
    return(FALSE)
  }
  # try for comma strings
  if(all(str_detect(v, "[0-9]+(,)[0-9]+"))) {
    v <- str_replace_all(v, ",", ".")
  }
  suppressWarnings(v_num <- as.numeric(v))
  na_orig <- sum(is.na(v))
  na_num <- sum(is.na(v_num))
  n_unique <- length(unique(v_num[!is.na(v_num)]))
  total <- length(v)
  (na_num - na_orig) / total < 0.01 && n_unique > 5
})]

if(length(numeric_like_columns) > 0) {
  print(numeric_like_columns)
  data_df <- data_df %>%
    mutate(across(all_of(numeric_like_columns), ~{
      # Check if column contains comma-separated numbers
      if(all(str_detect(na.omit(.), "[0-9]+(,)[0-9]+"))) {
        # Replace comma with period, then convert
        as.numeric(str_replace_all(., ",", "."))
      } else {
        # Direct conversion
        as.numeric(.)
      }
    }))
}

numeric_vars <- names(data_df)[sapply(data_df, is.numeric)]

## use distributions to identify outliers

# 1. For each numeric variable, compute IQR bounds and collect outlier row‐indices
outliers <- map_dfr(numeric_vars, function(var) {
  x <- data_df[[var]]
  # skip too‐sparse or constant vectors
  if (length(unique(na.omit(x))) < 10) {
    return(NULL)
  }
  qnt <- quantile(x, c(0.05, 0.95), na.rm = TRUE)
  iqr <- qnt[2] - qnt[1]
  lb <- qnt[1] - 1 * iqr
  ub <- qnt[2] + 1 * iqr
  idx <- which(x < lb | x > ub)
  if (length(idx) == 0) {
    return(NULL)
  }
  tibble(
    variable    = var,
    row         = idx,
    value       = x[idx],
    lower_bound = lb,
    upper_bound = ub
  )
})



# 2. Summary: how many outliers per variable
outliers %>%
  count(variable, name = "n_outliers") %>%
  mutate(pct = n_outliers / nrow(data_df) * 100) %>%
  arrange(desc(n_outliers)) %>%
  print()

# 3. Optional: get the set of row‐numbers to drop (any variable)
outlier_rows <- unique(outliers$row)

# 4. Plot each variable with outliers highlighted
output_dir <- "discrepancies/2025-10-Data-Version/outliers/"
dir.create(output_dir, showWarnings = FALSE)
for (var in unique(outliers$variable)) {
  df <- tibble(
    row = seq_len(nrow(data_df)),
    value = data_df[[var]]
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
      x = "Row index", y = "Value"
    ) +
    theme_bw()
  ggsave(filename = paste0(output_dir,"outliers_", gsub("[^A-Za-z0-9]", "_", var), ".png"), plot = p, width = 12, height = 4)
}

# order factors

quality_order <- c("Poor", "Satisfactory", "Good")
absent_present_order <- c("Absent", "Minimal", "Prominent")
no_yes_order <- c("No", "Yes")
lowhigh_order <- c("Low", "High")
lowhigh_order_2 <- c("No", "Low", "Intermediate", "High")
absent_present_order_2 <- c("Absent", "Minimal", "Moderate", "Extensive")

orders_list <- list(quality_order, absent_present_order, absent_present_order_2, no_yes_order, lowhigh_order, lowhigh_order_2)

for (var in names(data_df)[sapply(data_df, is.character)]) {
  vals <- na.omit(unique(data_df[[var]]))
  # try each ordering; once matched, convert and move on
  for (ord in orders_list) {
    if (length(vals) > 0 && all(vals %in% ord)) {
      data_df[[var]] <- factor(data_df[[var]], levels = ord, ordered = TRUE)
      break
    }
  }
}

unordered_vars <- names(data_df)[sapply(data_df, function(v) {
  (is.character(v) || is.factor(v)) && !is.ordered(v)
})]

if (length(unordered_vars) > 0) {
  message("Unordered character/factor variables:")
  print(unordered_vars)
} else {
  message("No unordered character/factor variables found")
}

## plots of distributions

# Create a directory for distribution plots if it doesn't exist
output_dir <- "discrepancies/2025-10-Data-Version/distributions/"
dir.create(output_dir, showWarnings = FALSE)

# Plot distributions for each variable in data_df
for (var in names(data_df)) {
  p <- ggplot(data_df, aes_string(x = paste0("`",var,"`"))) +
    theme_bw() +
    labs(title = paste("Distribution of", var))
  
  if (is.numeric(data_df[[var]])) {
    p <- p + geom_histogram(bins=50, fill = "blue", color = "black", alpha = 0.7)
  } else if (is.factor(data_df[[var]]) || is.character(data_df[[var]])) {
    p <- p + geom_bar(fill = "blue", color = "black", alpha = 0.7)
  } else {
    next  # Skip if the variable is neither numeric nor categorical
  }
  
  ggsave(filename = paste0(output_dir, "distribution_", gsub("[^A-Za-z0-9]", "_", var), ".png"), plot = p, width = 12, height = 6)
}
na_columns <- colnames(data_df)[sapply(data_df, function(x) any(is.na(x)))]
if (length(na_columns) > 0) {
  message("Columns with NAs found:")
  print(na_columns)
} else {
  message("No columns with NAs found.")
}

# Find and remove rows with NA values
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

data_df_pre_scaling <- data_df
dir.create("discrepancies/2025-10-Data-Version/processed_data/", showWarnings = FALSE)
saveRDS(data_df_pre_scaling, file = "discrepancies/2025-10-Data-Version/processed_data/data_df_pre_scaling.rds")
# scaling numeric variables
scale_ <- function(x) {
  return(scale(x, center = TRUE, scale=TRUE) %>% as.numeric())
}
data_df <- data_df %>%
  mutate(across(all_of(numeric_vars), scale_))

saveRDS(data_df, file = "discrepancies/2025-10-Data-Version/processed_data/data_df.rds")

variables <- variables %>% mutate(interaction="")
saveRDS(variables, file = "discrepancies/2025-10-Data-Version/processed_data/variables.rds")