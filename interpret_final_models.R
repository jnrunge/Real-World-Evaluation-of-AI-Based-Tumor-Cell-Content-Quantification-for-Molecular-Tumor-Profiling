## p value of model terms

# Helper to remove terms and compare models
compare_model_drop_terms <- function(model, drop_terms) {
    reduced_formula <- update(
        formula(model),
        as.formula(
            paste(". ~ . -", paste(drop_terms, collapse = " - "))
        )
    )
    reduced_model <- lm(reduced_formula, data = model$model)
    anova_res <- anova(reduced_model, model)
    list(
        reduced_formula = reduced_formula,
        reduced_model = reduced_model,
        anova = anova_res
    )
}