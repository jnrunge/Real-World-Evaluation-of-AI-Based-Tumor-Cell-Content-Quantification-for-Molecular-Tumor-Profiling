### would hemorraghe be btter than any of the other variables?

data <- all_model_variables_renamed

fit <- all_types_not_enforced_interaction
tt <- terms(fit)
term_labels <- attr(tt, "term.labels")
response_str <- deparse(formula(fit)[[2]])
intercept_flag <- as.integer(attr(tt, "intercept"))
new_var <- "Hemorrhage.Resection.Biopsy._Extensive.Blood.Content..Cytology."
base_AIC <- AIC(fit)

# Try to capture any offset terms (if present)
offset_idx <- attr(tt, "offset")
offset_terms <- character(0)
if (!is.null(offset_idx) && length(offset_idx) > 0) {
    vars <- attr(tt, "variables")
    offset_terms <- vapply(offset_idx, function(i) deparse(vars[[i + 1]]), character(1))
}

build_formula <- function(terms_vec) {
    f <- reformulate(terms_vec, response = response_str, intercept = intercept_flag)
    if (length(offset_terms)) {
        rhs_terms <- attr(terms(f), "term.labels")
        rhs <- paste(c(rhs_terms, paste0("offset(", offset_terms, ")")), collapse = " + ")
        as.formula(paste(response_str, "~", if (intercept_flag == 0) "-1 + " else "", rhs))
    } else {
        f
    }
}

results <- lapply(term_labels, function(drop_term) {
    # Remove the term, add the new var (ensure it appears exactly once)
    new_terms <- unique(c(setdiff(term_labels, drop_term), new_var))
    f_new <- build_formula(new_terms)
    m_new <- try(suppressWarnings(update(fit, formula = f_new)), silent = TRUE)

    if (inherits(m_new, "try-error")) {
        data.frame(
            term_removed = drop_term,
            AIC_original = base_AIC,
            AIC_new = NA_real_,
            delta_AIC = NA_real_,
            fit_ok = FALSE,
            stringsAsFactors = FALSE
        )
    } else {
        a_new <- AIC(m_new)
        data.frame(
            term_removed = drop_term,
            AIC_original = base_AIC,
            AIC_new = a_new,
            delta_AIC = a_new - base_AIC,
            fit_ok = TRUE,
            stringsAsFactors = FALSE
        )
    }
})

res <- do.call(rbind, results)
res <- res[order(res$delta_AIC), ]
res




# Leave-one-out analysis: remove one row at a time, add new_var, fit model, get AIC change
data <- all_model_variables_renamed
fit <- all_types_not_enforced_interaction$model
tt <- terms(fit)
term_labels <- attr(tt, "term.labels")
new_var <- "Hemorrhage.Resection.Biopsy._Extensive.Blood.Content..Cytology."

# Get the original formula and update it to include new_var
orig_formula <- formula(fit)
new_terms <- unique(c(term_labels, new_var))
f_new <- build_formula(new_terms)

# Number of rows
n <- nrow(data)

# Function to fit model on subset and get AIC
fit_and_aic <- function(subset_data, formula) {
    m <- try(suppressWarnings(lm(formula, data = subset_data)), silent = TRUE)
    if (inherits(m, "try-error")) {
        NA_real_
    } else {
        AIC(m)
    }
}

# Parallel setup (using parallel package)
library(parallel)
num_cores <- detectCores() - 1 # Leave one core free

# Results list
results <- vector("list", 3) # For k=1,2,3

for (k in 1:3) {
    # Generate all combinations of k rows to remove
    combos <- combn(n, k, simplify = FALSE)

    # Function to process each combination
    process_combo <- function(combo) {
        subset_data <- data[-combo, ]
        # Compute base AIC on subsetted data with original formula
        base_AIC <- fit_and_aic(subset_data, orig_formula)
        # Compute new AIC on subsetted data with new formula
        aic_new <- fit_and_aic(subset_data, f_new)
        delta_aic <- aic_new - base_AIC
        data.frame(
            k = k,
            rows_removed = paste(combo, collapse = ","),
            AIC_original = base_AIC,
            AIC_new = aic_new,
            delta_AIC = delta_aic,
            fit_ok = !is.na(aic_new),
            stringsAsFactors = FALSE
        )
    }

    # Parallelize over combinations
    results[[k]] <- mclapply(combos, process_combo, mc.cores = num_cores)
}

# Combine results
res <- do.call(rbind, unlist(results, recursive = FALSE))
res <- res[order(res$delta_AIC), ]
res
