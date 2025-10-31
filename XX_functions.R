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