
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
    var
}


perform_forward_selection <- function(v_table, model_data, title_prefix = "", force_interaction=TRUE, delta_aic=2, steps=2000, seed=123, reps=1,
model_pre_done=NULL) {

if (is.null(model_pre_done)) {
    v_table <- v_table %>%
        mutate(
            variable = make.names(variable, unique = TRUE),
            interaction = if_else(
                interaction != "",
                sapply(strsplit(interaction, ";", fixed = TRUE), function(parts) {
                    paste(make.names(parts, unique = FALSE), collapse = ";")
                }),
                interaction
            )
        )
    colnames(model_data) <- make.names(colnames(model_data), unique = TRUE)


    # Convert any character variables in v_table to factors in model_data
    char_vars <- intersect(v_table$variable, names(model_data))
    char_vars <- char_vars[sapply(model_data[char_vars], is.character)]
    if (length(char_vars) > 0) {
        model_data <- model_data %>%
            mutate(across(all_of(char_vars), as.factor))
    }


    # Forward stepwise regression
    if (reps == 1) {
        fit <- stepwise_aic(model_data, v_table,
            direction = "both", trace = TRUE, steps = steps, enforce_interactions = force_interaction, aic_delta = delta_aic, seed = seed
        )
    } else {
        # run 'reps' fits in parallel using different random seeds
        set.seed(seed)
        seeds <- sample.int(.Machine$integer.max, reps)
        ncores <- min(reps, parallel::detectCores())
        library(doParallel)
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        fits <- foreach(
            s = seeds,
            .packages = c("tidyverse", "stats")
        ) %dopar% {
            set.seed(s)
            # perform the stepwise selection and return the model object
            source("/mnt/S/People/JanN/PathAI_analyses/discrepancies/manual-stepwise.R")
            stepwise_aic(
                model_data, v_table,
                direction = "both",
                trace = TRUE,
                steps = steps,
                enforce_interactions = force_interaction,
                aic_delta = delta_aic,
                seed = s
            )
        }
        stopCluster(cl)
        # collect results as a list of model fits
    }

    # if multiple reps were run, extract their AICs, plot the distribution, and pick the best
    if (reps > 1 && exists("fits")) {
        # 1. compute AIC for each returned fit
        aic_vals <- sapply(fits, function(xx) AIC(xx[[1]]))

        # 2. plot the distribution of AICs
        library(ggplot2)
        p_aic <- ggplot(data.frame(AIC = aic_vals), aes(x = AIC)) +
            geom_histogram(
                binwidth = diff(range(aic_vals)) / 30,
                fill = "steelblue", color = "white"
            ) +
            labs(
                title = "Distribution of AIC across fits",
                x = "AIC", y = "Count"
            ) +
            theme_minimal()
        print(p_aic)
        # 3. select the fit with the lowest AIC
        best_idx <- which.min(aic_vals)[1]
        fit <- fits[[best_idx]]
    }


    forward_model <- fit[[1]]
    steps <- fit[[2]]

    # Manually force‐include Sample.type in the final model
    # if (!"Sample.type" %in% names(coef(forward_model)) &&  "Sample type" %in% v_table$variable) {
    #    forward_model <- update(forward_model, . ~ . + Sample.type, data = model_data)
    # }

    plot(AIC ~ step, steps)
}

if(!is.null(model_pre_done)) {
    forward_model <- model_pre_done
}
    
set_theme(
    geom.outline.color    = "antiquewhite4",
    geom.outline.size     = 1,
    geom.label.size       = 2,
    geom.label.color      = "grey50",
    title.color           = "red",
    title.size            = 1,
    geom.errorbar.size    = 1.5,  # increase error bar thickness
    base = theme_bw()
)

    # predicted‐value plots
    pred_plots <- plot_model(
        forward_model,
        type      = "pred",
        show.data = TRUE,
        jitter    = 0.1,
    )

    # variable importance usign permuation importance

    var_importance_df<-vip(forward_model, num_features=1000)$data

    # Build a mapping from each model term back to its v_table$variable
    term_map <- tibble(
        model_term = var_importance_df$Variable
    ) %>%
    rowwise() %>%
    mutate(
        orig_vars = list({
            # Split interaction terms
            parts <- strsplit(model_term, ":", fixed = TRUE)[[1]]
            # For each part, find matching v_table$variable
            hits <- sapply(parts, function(part) {
                vhits <- v_table$variable[
                    str_detect(
                        part,
                        stringr::regex(
                            paste0("^", v_table$variable),
                            ignore_case = FALSE
                        )
                    )
                ]
                if (length(vhits) == 1) {
                    vhits[[1]]
                } else if (length(vhits) > 1) {
                    # Filter vhits to keep only those with the maximum nchar
                    maxlen <- max(nchar(vhits))
                    vhits_max <- vhits[nchar(vhits) == maxlen]
                    if (length(vhits_max) == 1) {
                        vhits_max[[1]]
                    } else {
                        print(part)
                        print(vhits_max)
                        stop("Too many matches in vip df for interaction part after filtering by max nchar")
                    }
                } else {
                    part
                }
            }, USE.NAMES = FALSE)
            hits
        })
    ) %>%
    ungroup()

    term_map_long <- term_map %>%
        unnest_longer(orig_vars)

    # Sum importance by original variable
    importance_summary <- var_importance_df %>%
        full_join(term_map_long, by=c("Variable"="model_term")) %>%
        group_by(orig_vars) %>%
        summarise(TotalImportance = sum(Importance), .groups = "drop") %>%
        mutate(RelativeImportance = TotalImportance / sum(TotalImportance)) %>%
        arrange(desc(RelativeImportance))

    # optional: print or return
    print(importance_summary)

    # Extract residuals
    res <- resid(forward_model)


    # Compute the model prediction when all predictors are “held constant”
    newdata_ref <- model_data %>%
        summarise(
            across(where(is.numeric), ~ 0),
            across(where(is.factor), ~ {
                lvls <- levels(.x)
                # Find the first level that actually exists in the data
                for (l in lvls) {
                    if (any(.x == l, na.rm = TRUE)) return(l)
                }
                # fallback: just use the first level
                lvls[1]
            })
        )

    pred_const <- predict(forward_model, newdata = newdata_ref)
    message("Predicted TC_Path_minus_TC_AI at constant settings: ", pred_const)


    for ( j in 1:length(pred_plots)){
# Get fitted values of just the effect of x1
# Prepare data and extract the first predictor name
x1_var <- names(pred_plots)[j]

# Compute the contribution of x1_var to the linear predictor,
# whether it is numeric or a factor (with multiple dummy columns).
mm <- model.matrix(formula(forward_model), data = model_data)
# find all columns in the design matrix that correspond to x1_var
cols <- grep(paste0("^", x1_var), colnames(mm))
# multiply those columns by their coefficients and sum
x1_term <- as.numeric(mm[, cols, drop = FALSE] %*% coef(forward_model)[cols])

# Get all model terms except those related to x1_var
# (creates properly aligned partial residuals with prediction line)
all_terms <- predict(forward_model, type = "terms")
other_vars_terms <- all_terms[, -grep(paste0("^", x1_var), colnames(all_terms)), drop = FALSE]
other_vars_effect <- rowSums(other_vars_terms)
intercept <- coef(forward_model)["(Intercept)"]

# Compute adjusted partial residuals
partial_resid_x1 <- fitted(forward_model) - other_vars_effect + res

# Add to data
model_data$partial_resid_x1 <- partial_resid_x1

# add partial‐residual points to the existing pred_plots
# Use "K" formatting for thousands on x axis if numeric
pred_plots[[x1_var]] <-
    pred_plots[[x1_var]] +
    geom_jitter(
        data = model_data,
        aes_string(x = x1_var, y = "partial_resid_x1"),
        inherit.aes = FALSE,
        color = "gray40",
        alpha = 1, shape = 1, width = 0.1, height = 0.1
    ) +
    ggtitle(
        paste0(
            get_pretty_name(x1_var),
            " [",
            round(importance_summary$RelativeImportance[importance_summary$orig_vars == x1_var] * 100, 1),
            "% RI]"
        )
    ) +
    xlab(NULL) +
    ylab(NULL) +
    (if (is.numeric(model_data[[x1_var]])) {
        scale_x_continuous(labels = function(x) {
            sapply(x, function(val) {
                if (is.na(val)) {
                    return(NA_character_)
                }
                if(abs(val)>=1e9) {
                    paste0(formatC(val / 1e9, format = "f", digits = 1), "B")
                } else if (abs(val) >= 1e6) {
                    paste0(formatC(val / 1e6, format = "f", digits = 1), "M")
                } else if (abs(val) >= 1e3) {
                    paste0(formatC(val / 1e3, format = "f", digits = ifelse(abs(val) >= 1e4, 0, 1)), "K")
                } else {
                    as.character(val)
                }
            })
        })
    } else {
        scale_x_discrete()
    })
    }
    

    # only attempt interaction plots if the model has “:” terms
    coef_nm <- names(coef(forward_model))
    if (any(grepl(":", coef_nm))) {
        warning("Interactions still not fully supported in plots optimization, but will be shown anyway.")
        # reorder interaction terms in the model so that any numeric variable comes first
        # for int plots with "many levels" because numeric I tested this fix
#         reorder_interaction_terms <- function(model, model_data) {
#             # Extract original formula and its terms
#             orig_f  <- formula(model)
#             terms_f <- attr(terms(orig_f), "term.labels")

#             # Reorder each interaction so numeric vars come first
#             new_terms <- vapply(terms_f, function(t) {
#                 if (!grepl(":", t, fixed = TRUE)) return(t)
#                 parts    <- strsplit(t, ":", fixed = TRUE)[[1]]
#                 is_num   <- vapply(parts, function(v) is.numeric(model_data[[v]]), logical(1))
#                 parts    <- parts[order(!is_num)]
#                 paste(parts, collapse = ":")
#             }, character(1))

#             # Build a new formula and refit
#             lhs <- deparse(orig_f[[2]])
#             rhs <- paste(new_terms, collapse = " + ")
#             new_f <- as.formula(paste(lhs, "~", rhs))
#             lm(new_f, data = model_data)
#         }

#         # Usage inside your chunk:
#         forward_model <- reorder_interaction_terms(forward_model, model_data)

#         # randomly drop one interaction term to create a second model
#         # set how many interaction terms you want to keep
#         # Repeat the random dropping and plotting 100 times, record errors and interaction terms

#         results <- tibble(
#             iteration   = integer(),
#             inter_terms = character(),
#             error       = logical(),
#             error_msg   = character()
#         )
# set.seed(123)
# library(foreach)
# library(doParallel)
# # use up to 16 cores
# n_cores <- min(16, parallel::detectCores())
# cl <- makeCluster(n_cores)
# registerDoParallel(cl)

# # parallel loop: each iteration returns a tibble row
# results <- foreach(
#     iter = seq_len(100),
#     .combine = bind_rows,
#     .packages = c("sjPlot", "dplyr")
# ) %dopar% {
#     # identify interaction terms in the current model
#     all_terms <- attr(terms(forward_model), "term.labels")
#     inter_terms <- grep(":", all_terms, value = TRUE)

#     # randomly keep all but one interaction term
#     n_keep <- length(inter_terms)-1
#     keep_terms <- if (n_keep > 0) sample(inter_terms, n_keep) else character(0)
#     drop_terms <- setdiff(inter_terms, keep_terms)

#     # update model by dropping the other interactions
#     forward_model2 <- if (length(drop_terms) > 0) {
#         frm <- as.formula(paste(". ~ . -", paste(drop_terms, collapse = " - ")))
#         reorder_interaction_terms(update(forward_model, frm, data = model_data), model_data)
#     } else {
#         forward_model
#     }

#     # attempt to get interaction plots, catch errors
#     err_flag <- FALSE
#     err_msg <- NA_character_
#     int_plots_tmp <- tryCatch(
#         {
#             plot_model(forward_model2, type = "int")
#         },
#         error = function(e) {
#             err_flag <<- TRUE
#             err_msg <<- e$message
#             NULL
#         }
#     )

#     # return one-row tibble
#     tibble(
#         iteration   = iter,
#         inter_terms = paste(attr(terms(forward_model2), "term.labels") %>% grep(pattern=":", value = TRUE), collapse = ","),
#         error       = err_flag,
#         error_msg   = err_msg
#     )
# }

# stopCluster(cl)

#         # Print the summary of which iterations failed and with which terms
#         results %>%
#             separate_rows(inter_terms, sep = ",") %>%
#             group_by(inter_terms) %>%
#             summarize(e_rate = sum(error) / n()) %>%
#             arrange(desc(e_rate))


# Generate interaction effect plots
all_terms <- attr(terms(forward_model), "term.labels")
int_terms <- grep(":", all_terms, value = TRUE)
int_plots <- list()
if (length(int_terms) > 0) {
    for (tm in int_terms) {
        parts <- strsplit(tm, ":", fixed = TRUE)[[1]]
        # detect which parts are numeric in model_data
        is_num <- sapply(parts, function(x) is.numeric(model_data[[x]]))
        # ensure the numeric term comes first
        parts <- parts[order(!is_num)]
        if (all(is_num)) {
            # both numeric: pick 4 quantile cut-points for the second variable
            sec <- parts[2]
            qs <- quantile(model_data[[sec]], probs = seq(0, 1, length.out = 4), na.rm = TRUE)
            sec_term <- sprintf("%s[%s]", sec, paste(round(qs, 2), collapse = ","))
            terms_arg <- c(parts[1], sec_term)
        } else {
            # at least one factor: just feed the two terms in order
            terms_arg <- parts
        }
        int_plots[[tm]] <- plot_model(
            forward_model,
            type  = "pred",
            terms = terms_arg,
            show.data = TRUE,
            jitter = 0.1
        )
    }
}
    

#         if(length(int_plots)==11 && names(int_plots)[11]=="labels"){
#             # dann ist length eigentlich 1 und die logik des durch die namen gehens geht nicht

#             # get the name of the interaction vars
#             formula(forward_model) %>% as.character() -> explanatory_variables_from_formula
#             terms_forward_model <- explanatory_variables_from_formula[3] %>% strsplit(split = " + ", fixed = TRUE)
#             terms_forward_model <- unlist(terms_forward_model)
#             interaction_var<-terms_forward_model[grepl(":", terms_forward_model, fixed = TRUE)]
            
#             int_plots <- list(int1 = int_plots)
#             names(int_plots) <- interaction_var
#         }else{
#             # otherwise, we have a list of interaction plots

#             # inference of interaction variable names

#             for(ints in seq_len(length(int_plots))){
#                 v1 <- int_plots[[ints]]$labels$x
#                 v2 <- int_plots[[ints]]$labels$colour

# formula(forward_model) %>% as.character() -> explanatory_variables_from_formula
# terms_forward_model <- explanatory_variables_from_formula[3] %>% strsplit(split = " + ", fixed = TRUE)
# terms_forward_model <- unlist(terms_forward_model)
# interaction_vars <- terms_forward_model[grepl(":", terms_forward_model, fixed = TRUE)]

# # build the two possible interaction strings
# target1 <- paste(v1, v2, sep = ":")
# target2 <- paste(v2, v1, sep = ":")

# # compute minimal edit‐distance of each candidate in interaction_vars to either target
# dists <- sapply(interaction_vars, function(iv) {
#     min(adist(iv, target1), adist(iv, target2))
# })

# # pick the best match
# best_match <- interaction_vars[which.min(dists)]

# # assign that name to the current interaction plot
# names(int_plots)[ints] <- best_match

#             }
#         }


    } else {
        int_plots <- list()
    }

    # ensure each is a list of ggplots and then concatenate
    pred_list <- if (inherits(pred_plots, "list")) pred_plots else list(pred_plots)
    int_list  <- if (inherits(int_plots,  "list")) int_plots  else list(int_plots)

    partial_plots_unique <- c(pred_list, int_list)


    

    if (!is.null(model_pre_done)) {
        return(list(
            plots = partial_plots_unique,
            importance_summary = importance_summary
        ))
    }else{
return(list(
    model = forward_model,
    plots = partial_plots_unique,
    steps = steps,
    importance_summary = importance_summary,
    fits = if (exists("fits")) fits else NULL,
    best_idx = if (exists("best_idx")) best_idx else NULL
))
    }
}
