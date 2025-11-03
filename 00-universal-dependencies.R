library(tidyverse)

pretty_names <- read_csv("discrepancies/variable_pretty_names.csv")
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
    write_csv(tibble(variable = var, pretty_name = var), "discrepancies/variable_pretty_names.csv", append = TRUE)
    var
}

response_vars <- variables %>%
    filter(type == "response") %>%
    pull(variable)
