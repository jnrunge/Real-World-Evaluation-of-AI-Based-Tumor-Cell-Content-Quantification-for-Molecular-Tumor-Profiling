# Optional dependency: MASS is not used by the function below (base stats::AIC/lm are used).
# Keep this only if you plan to use MASS::stepAIC elsewhere in this script.
library(MASS) # 7.3-60.0.1

#' Randomized stepwise AIC model selection with optional enforced interactions
#'
#' This function performs a randomized add/remove search over candidate terms
#' (main effects and optionally interactions) using AIC as the criterion.
#' It can "enforce" interactions so that a variable, its interacting partner,
#' and their interaction enter/leave the model together.
#'
#' Note:
#' - direction and trace are accepted for interface compatibility, but currently ignored.
#' - Rows with missing values in required variables are dropped before fitting.
#' - MASS is not used here; only base stats::lm and stats::AIC are used.
#'
#' @param model_data A data.frame with all candidate variables as columns.
#' @param v_table    A data.frame with columns:
#'                   - variable: names matching colnames(model_data)
#'                   - type:     "response" or "explanatory"
#'                   - interaction: "" or a semicolon-separated list of variable
#'                     names that this variable must/can interact with
#' @param direction  Ignored; reserved for future compatibility ("both"|"forward"|"backward").
#' @param trace      Ignored; reserved for future compatibility.
#' @param steps      Integer; maximum number of randomized add/remove proposals.
#' @param aic_delta  Numeric; minimal AIC improvement required to accept a change.
#' @param enforce_interactions Logical; if TRUE, a variable, its interacting partner(s),
#'                   and their interaction(s) are grouped and added/removed together.
#'                   If FALSE, each term (main effect and interaction) is considered independently.
#' @param seed       Integer; RNG seed for reproducibility of the randomized search.
#' @return           A list with:
#'                     - model: final stats::lm object after selection
#'                     - steps: data.frame logging the search (see note on column names below)
#' @details
#' - Interactions in v_table$interaction may be provided as "varA;varB" to indicate multiple partners.
#' - Required packages: if enforce_interactions = FALSE and removal branch is reached,
#'   code paths use stringr/magrittr helpers (str_replace, str_replace_all, %>%).
#'   Load those packages or adapt the code if you hit those paths.
#' - Known caveats retained for reproducibility:
#'   (1) Add step accepts only if AIC improves by > aic_delta AND the proposed AIC > 0,
#'       which rejects improvements to negative AIC; adjust if needed.
#'   (2) steps log is initialized with column 'terms' but rows use 'term'; align if you rely on it.
#'   (3) When nrow(raw_interactions) == 0, interaction_rows <- tibble() requires the tibble package;
#'       use data.frame() instead if tibble is unavailable.
stepwise_aic <- function(model_data, v_table, direction = "both", trace = FALSE, steps=1000, aic_delta=2, enforce_interactions=TRUE, seed=123) {
  # Filter out variables with <= 1 unique non-missing value (uninformative predictors)
  vars <- v_table$variable
  uniq_counts <- sapply(model_data[vars], function(x) length(unique(x[!is.na(x)])))
  bad_vars <- names(uniq_counts)[uniq_counts <= 1]
  if (length(bad_vars)) {
    v_table <- v_table[!v_table$variable %in% bad_vars, ]
  }

  # Identify response (must be exactly one) and main explanatory variables
  resp_vars <- v_table$variable[v_table$type == "response"]
  if (length(resp_vars) != 1) {
    stop("var_table must contain exactly one response variable.")
  }
  response <- resp_vars
  expl_vars <- v_table$variable[v_table$type == "explanatory"]

  # Build interaction rows, allowing multiple partners per variable via ";" in v_table$interaction
  raw_interactions <- v_table[v_table$type == "explanatory" & v_table$interaction != "", ]
  if (nrow(raw_interactions) > 0) {
    interaction_rows <- do.call(rbind, lapply(seq_len(nrow(raw_interactions)), function(i) {
      var  <- raw_interactions$variable[i]
      ints <- strsplit(raw_interactions$interaction[i], ";", fixed = TRUE)[[1]]
      data.frame(
        variable    = var,
        interaction = trimws(ints),
        stringsAsFactors = FALSE
      )
    }))
    # For display; actual quoted terms are constructed below
    interaction_terms <- paste0(interaction_rows$variable, ":", interaction_rows$interaction)
  } else {
    # Requires the tibble package; replace with data.frame() if tibble is unavailable in your environment.
    interaction_rows <- tibble()
    interaction_terms <- character(0)
  }

  # Preprocess: drop rows with NA in any required variable (response, all main effects, and any interacting partners)
  req_vars <- unique(c(
    response,
    expl_vars,
    if (nrow(interaction_rows) > 0) interaction_rows$interaction else character()
  ))
  orig_n   <- nrow(model_data)
  data     <- model_data[stats::complete.cases(model_data[, req_vars]), ]
  n_removed <- orig_n - nrow(data)  # kept for diagnostics; not returned

  # Assemble formula components with quoted names to handle spaces/special characters
  quote_name <- function(x) paste0("`", x, "`")
  response_q <- quote_name(response)
  main_terms_q <- sapply(expl_vars, quote_name)
  interaction_terms_q <- if (length(interaction_terms) > 0) {
    # Quote both sides of each interaction
    apply(interaction_rows, 1, function(r) {
      paste0(quote_name(r["variable"]), ":", quote_name(r["interaction"]))
    })
  } else {
    character(0)
  }
  full_terms_q  <- c(main_terms_q, interaction_terms_q)
  full_formula_str  <- paste(response_q, "~", paste(full_terms_q, collapse = " + "))
  empty_formula_str <- paste(response_q, "~ 1")

  full_formula  <- as.formula(full_formula_str)
  empty_formula <- as.formula(empty_formula_str)

  # Fit empty and full models on the cleaned data. lm_full is not used later, but can help validate term availability.
  lm_empty <- lm(empty_formula, data = data)
  lm_full  <- lm(full_formula,  data = data)

  # Define candidate "groups" for add/remove:
  # - If enforce_interactions = TRUE and a variable v has an interacting partner inter,
  #   we group: v (main), inter (main), and v:inter (interaction).
  # - If enforce_interactions = FALSE, main effects and interactions are considered independently.
  #   When testing interaction addition, we add only the minimal required terms (both mains if absent, or missing main + interaction).
  if (enforce_interactions) {
    # Original grouped behavior: everything moves together
    groups <- lapply(expl_vars, function(v) {
      if (v %in% interaction_rows$variable) {
        inter <- interaction_rows$interaction[interaction_rows$variable == v]
        c(paste0("`", v, "`"), paste0("`", v, "` * `", inter, "` + `", inter, "`"))
      } else {
        paste0("`", v, "`")
      }
    })
    names(groups) <- expl_vars
  } else {
    # Independent testing: build separate groups for main effects and interactions
    # Main effects groups
    groups <- lapply(expl_vars, function(v) {
      paste0("`", v, "`")
    })
    names(groups) <- expl_vars
    
    # Add interaction groups (these will be tested separately)
    if (nrow(interaction_rows) > 0) {
      for (i in seq_len(nrow(interaction_rows))) {
        v <- interaction_rows$variable[i]
        inter <- interaction_rows$interaction[i]
        int_name <- paste0(v, ":", inter)
        # Store both the interaction term and what main effects it requires
        groups[[int_name]] <- list(
          main_v = paste0("`", v, "`"),
          main_inter = paste0("`", inter, "`"),
          interaction = paste0("`", v, "`:`", inter, "`")
        )
      }
    }
  }

  # Initialize state and RNG
  current_terms <- character(0)
  best_aic      <- AIC(lm(empty_formula, data = data))
  set.seed(seed)

  # Steps log. Note: initialized with column 'terms', but rows below bind a 'term' column.
  # Harmonize if you plan to parse this table programmatically.
  steps_tbl <- data.frame(
    step = integer(0),
    terms = I(list()),
    AIC = numeric(0),
    Note = character(0),
    stringsAsFactors = FALSE
  )

  # Randomized add/remove loop
  for (i in seq_len(steps)) {
    rand_val <- runif(1)
    
    if (rand_val < 0.33 && length(current_terms) < length(groups)*2) {
      # ADDITION STEP
      # Determine which groups can be added (not already fully represented)
      if (enforce_interactions) {
        avail <- setdiff(names(groups), sub(":.*", "", current_terms))
      } else {
        # For independent mode, check each group separately
        avail <- character(0)
        for (g_name in names(groups)) {
          g <- groups[[g_name]]
          if (is.list(g)) {
            # This is an interaction group - check if interaction term is already present
            if (!(g$interaction %in% current_terms)) {
              avail <- c(avail, g_name)
            }
          } else {
            # This is a main effect - check if it's already present
            if (!(g %in% current_terms)) {
              avail <- c(avail, g_name)
            }
          }
        }
      }
      
      if (length(avail)) {
        g_name <- sample(avail, 1)
        g <- groups[[g_name]]
        
        # Determine what terms to add
        if (enforce_interactions || !is.list(g)) {
          # Simple case: add the group as-is
          terms_to_add <- g
          cand <- c(current_terms, terms_to_add)
        } else {
          # Interaction in independent mode: add only what's missing
          # Check which main effects are already in the model
          has_main_v <- g$main_v %in% current_terms
          has_main_inter <- g$main_inter %in% current_terms
          
          terms_to_add <- character(0)
          if (!has_main_v) {
            terms_to_add <- c(terms_to_add, g$main_v)
          }
          if (!has_main_inter) {
            terms_to_add <- c(terms_to_add, g$main_inter)
          }
          terms_to_add <- c(terms_to_add, g$interaction)
          
          cand <- c(current_terms, terms_to_add)
        }
        
        # Fit candidate model and test AIC improvement
        fm_str <- paste(response_q, "~", if (length(cand)) paste(cand, collapse=" + ") else "1")
        m <- lm(as.formula(fm_str), data = data)
        a <- AIC(m)
        # Accept only if AIC improves by > aic_delta AND proposed AIC > 0.
        # The "> 0" prevents accepting negative AIC even when better; adjust if you expect negative AIC.
        if ((distance <- (best_aic - a)) > aic_delta && a > 0) {
          current_terms <- cand
          best_aic <- a
          steps_tbl <- rbind(steps_tbl, data.frame(
            step = i,
            term = g_name,
            AIC  = best_aic,
            Note = paste0("add (Δ=", round(distance, 2), ", added: ", paste(terms_to_add, collapse = ", "), ")")
          ))
        } else {
          steps_tbl <- rbind(steps_tbl, data.frame(
            step = i,
            term = g_name,
            AIC  = best_aic,
            Note = paste0("skip add (candidate AIC=", round(a, 2), ", Δ=", round(best_aic - a, 2), ")")
          ))
        }
      }
    } else if (rand_val < 0.66 && length(current_terms)) {
      # REMOVAL STEP
      if (enforce_interactions) {
        # Original grouped removal behavior
        pres <- unique(sub(":.*", "", current_terms))
        pres <- pres[gsub("`", "", pres) %in% names(groups)]
        if (length(pres) == 0) next
        
        g_name <- sample(pres, 1)
        g_name_clean <- gsub("`", "", g_name)
        terms_to_remove <- groups[[g_name_clean]]
        cand <- setdiff(current_terms, terms_to_remove)
      } else {
        # Independent removal: identify what can be removed
        removable <- character(0)
        for (g_name in names(groups)) {
          g <- groups[[g_name]]
          if (is.list(g)) {
            # Interaction group: can remove if interaction is present
            if (g$interaction %in% current_terms) {
              removable <- c(removable, g_name)
            }
          } else {
            # Main effect: can remove if present AND not required by any interaction in the model
            if (g %in% current_terms) {
              # Check if this main effect is part of any interaction currently in model
              is_required <- FALSE
              for (other_name in names(groups)) {
                other <- groups[[other_name]]
                if (is.list(other) && other$interaction %in% current_terms) {
                  if (other$main_v == g || other$main_inter == g) {
                    is_required <- TRUE
                    break
                  }
                }
              }
              if (!is_required) {
                removable <- c(removable, g_name)
              }
            }
          }
        }
        
        if (length(removable) == 0) next
        
        g_name <- sample(removable, 1)
        g <- groups[[g_name]]
        
        # Determine what to remove
        if (is.list(g)) {
          # Removing interaction: remove only the interaction term itself
          terms_to_remove <- g$interaction
        } else {
          # Removing main effect
          terms_to_remove <- g
        }
        
        cand <- setdiff(current_terms, terms_to_remove)
      }
      
      # Fit candidate model and test if removal is acceptable
      fm_str <- paste(response_q, "~", if (length(cand)) paste(cand, collapse=" + ") else "1")
      m <- lm(as.formula(fm_str), data = data)
      a <- AIC(m)
      
      # Accept removal if AIC does not worsen by >= aic_delta
      if ((distance <- (a - best_aic)) < aic_delta) {
        current_terms <- cand
        best_aic <- a
        steps_tbl <- rbind(steps_tbl, data.frame(
          step = i,
          term = g_name,
          AIC  = best_aic,
          Note = paste0("remove (Δ=", round(distance, 2), ", removed: ", paste(terms_to_remove, collapse=", "), ")")
        ))
      } else {
        steps_tbl <- rbind(steps_tbl, data.frame(
          step = i,
          term = g_name,
          AIC  = best_aic,
          Note = paste0("keep (removal AIC=", round(a, 2), ", Δ=", round(a - best_aic, 2), ")")
        ))
      }
    } else if (length(current_terms) > 0) {
      # REPLACEMENT STEP
      # Identify what can be removed (same logic as removal)
      if (enforce_interactions) {
        pres <- unique(sub(":.*", "", current_terms))
        pres <- pres[gsub("`", "", pres) %in% names(groups)]
      } else {
        # Independent removal: identify what can be removed
        pres <- character(0)
        for (g_name in names(groups)) {
          g <- groups[[g_name]]
          if (is.list(g)) {
            # Interaction group: can remove if interaction is present
            if (g$interaction %in% current_terms) {
              pres <- c(pres, g_name)
            }
          } else {
            # Main effect: can remove if present AND not required by any interaction in the model
            if (g %in% current_terms) {
              # Check if this main effect is part of any interaction currently in model
              is_required <- FALSE
              for (other_name in names(groups)) {
                other <- groups[[other_name]]
                if (is.list(other) && other$interaction %in% current_terms) {
                  if (other$main_v == g || other$main_inter == g) {
                    is_required <- TRUE
                    break
                  }
                }
              }
              if (!is_required) {
                pres <- c(pres, g_name)
              }
            }
          }
        }
      }
      
      # Identify what can be added (same logic as addition)
      if (enforce_interactions) {
        avail <- setdiff(names(groups), sub(":.*", "", current_terms))
      } else {
        avail <- character(0)
        for (g_name in names(groups)) {
          g <- groups[[g_name]]
          if (is.list(g)) {
            if (!(g$interaction %in% current_terms)) {
              avail <- c(avail, g_name)
            }
          } else {
            if (!(g %in% current_terms)) {
              avail <- c(avail, g_name)
            }
          }
        }
      }
      
      # Proceed only if we have both something to remove and something to add
      if (length(pres) > 0 && length(avail) > 0) {
        # Sample one term to remove and one to add
        g_name_remove <- sample(pres, 1)
        g_name_add <- sample(avail, 1)
        
        # Determine terms to remove
        if (enforce_interactions) {
          g_name_clean <- gsub("`", "", g_name_remove)
          terms_to_remove <- groups[[g_name_clean]]
        } else {
          g_remove <- groups[[g_name_remove]]
          if (is.list(g_remove)) {
            terms_to_remove <- g_remove$interaction
          } else {
            terms_to_remove <- g_remove
          }
        }
        
        # Remove the old terms first
        cand_after_remove <- setdiff(current_terms, terms_to_remove)
        
        # Determine terms to add
        g_add <- groups[[g_name_add]]
        if (enforce_interactions || !is.list(g_add)) {
          terms_to_add <- g_add
          cand <- c(cand_after_remove, terms_to_add)
        } else {
          # Interaction in independent mode: add only what's missing
          has_main_v <- g_add$main_v %in% cand_after_remove
          has_main_inter <- g_add$main_inter %in% cand_after_remove
          
          terms_to_add <- character(0)
          if (!has_main_v) {
            terms_to_add <- c(terms_to_add, g_add$main_v)
          }
          if (!has_main_inter) {
            terms_to_add <- c(terms_to_add, g_add$main_inter)
          }
          terms_to_add <- c(terms_to_add, g_add$interaction)
          
          cand <- c(cand_after_remove, terms_to_add)
        }
        
        # Fit candidate model and test AIC improvement
        fm_str <- paste(response_q, "~", if (length(cand)) paste(cand, collapse=" + ") else "1")
        m <- lm(as.formula(fm_str), data = data)
        a <- AIC(m)
        
        # Accept replacement if AIC improves by > aic_delta AND proposed AIC > 0
        if ((distance <- (best_aic - a)) > aic_delta && a > 0) {
          current_terms <- cand
          best_aic <- a
          steps_tbl <- rbind(steps_tbl, data.frame(
            step = i,
            term = paste0(g_name_remove, "→", g_name_add),
            AIC  = best_aic,
            Note = paste0("replace (Δ=", round(distance, 2), ", removed: ", 
                         paste(terms_to_remove, collapse=", "), ", added: ", 
                         paste(terms_to_add, collapse=", "), ")")
          ))
        } else {
          steps_tbl <- rbind(steps_tbl, data.frame(
            step = i,
            term = paste0(g_name_remove, "→", g_name_add),
            AIC  = best_aic,
            Note = paste0("skip replace (candidate AIC=", round(a, 2), ", Δ=", 
                         round(best_aic - a, 2), ", would remove: ", 
                         paste(terms_to_remove, collapse=", "), ", would add: ", 
                         paste(terms_to_add, collapse=", "), ")")
          ))
        }
      }
    }
  }

  # Build final model from the selected terms (whitespace around '+' is inconsequential for parsing).
  final_model <- lm(
    as.formula(paste(response_q, "~",
      if (length(current_terms)) paste(current_terms, collapse=" +") else "1"
    )),
    data = data
  )

  # Return the fitted model and the logged steps.
  return(list(
    model = final_model,
    steps = steps_tbl
  ))
}

# Example usage (commented):
# - Ensure any needed packages (stringr, tibble) are available if you hit those code paths.
# - Provide v_table with columns variable, type, interaction (semicolon-separated for multiple).
