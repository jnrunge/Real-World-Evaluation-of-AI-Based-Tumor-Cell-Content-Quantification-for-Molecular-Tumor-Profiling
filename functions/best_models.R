perform_forward_selection <- function(
    v_table, model_data, title_prefix = "", force_interaction = TRUE, delta_aic = 2, steps = 2000, seed = 123, reps = 1,
    model_pre_done = NULL) {
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
            cl <- makeCluster(ncores)
            registerDoParallel(cl)
            fits <- foreach(
                s = seeds,
                .packages = c("tidyverse", "stats"),
                .export = "project_dir"
            ) %dopar% {
                set.seed(s)
                # perform the stepwise selection and return the model object
                source(file.path(project_dir, "functions/manual-stepwise.R"))
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

        #     if (reps > 1 && exists("fits")) {
        #         # Combine all steps data from fits[[]][[2]] into one tibble
        #         steps_combined <- bind_rows(
        #             lapply(seq_along(fits), function(i) {
        #                 fits[[i]][[2]] %>%
        #                     mutate(fit_index = i)
        #             })
        #         )

        #         # Plot AIC ~ step faceted by fit_index
        #         p_steps <- ggplot(steps_combined, aes(x = step, y = AIC)) +
        #             geom_line() +
        #             geom_point(size = 0.5) +
        #             facet_wrap(~ fit_index, scales = "free_y") +
        #             labs(
        #                 title = "AIC progression across stepwise selection",
        #                 x = "Step",
        #                 y = "AIC"
        #             ) +
        #             theme_bw()()

        #     } else {
        #         p_steps <- plot(AIC ~ step, steps)
        #     }
        #
    }

    if (!is.null(model_pre_done)) {
        forward_model <- model_pre_done
    }


    # variable importance usign permuation importance

    var_importance_df <- vip(forward_model, num_features = 1000)$data

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
        full_join(term_map_long, by = c("Variable" = "model_term")) %>%
        group_by(orig_vars) %>%
        summarise(TotalImportance = sum(Importance), .groups = "drop") %>%
        mutate(RelativeImportance = TotalImportance / sum(TotalImportance)) %>%
        arrange(desc(RelativeImportance))

    




    if (!is.null(model_pre_done)) {
        return(list(
            importance_summary = importance_summary
        ))
    } else {
        return(list(
            model = forward_model,
            steps = steps,
            importance_summary = importance_summary,
            fits = if (exists("fits")) fits else NULL,
            best_idx = if (exists("best_idx")) best_idx else NULL
        ))
    }
}
