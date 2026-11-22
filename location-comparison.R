# Run the setup, load-data, hours-to-close, and cleaning chunks in
# "311-cases - Caterina.Rmd" first. Then Source this file in RStudio.
# This script uses base R and prints the model comparison in the console.

compare_location <- function(cases_clean, rare_cutoff = 25L) {
  location <- c("police_district", "supervisor_district", "analysis_neighborhood")
  baseline <- c("source", "service_name")
  required <- c("opened", "closed", "hours_to_close", baseline, location)
  missing_columns <- setdiff(required, names(cases_clean))
  if (length(missing_columns)) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }

  data <- cases_clean[, required]
  data <- data[!is.na(data$opened) & !is.na(data$closed) &
                 is.finite(data$hours_to_close) & data$hours_to_close >= 0, ]
  data <- data[order(data$opened), ]
  opening_times <- unique(data$opened)
  if (length(opening_times) < 3L) {
    stop("Need at least three distinct opening timestamps.")
  }

  # Use roughly the first 60%, next 20%, and last 20% by opening time.
  # All cases with the same opening timestamp stay on the same side.
  last <- nrow(data)
  validation_start <- min(data$opened[floor(0.60 * last) + 1L],
                          opening_times[length(opening_times) - 1L])
  validation_start <- max(validation_start, opening_times[2L])
  test_start <- max(data$opened[floor(0.80 * last) + 1L],
                    opening_times[which(opening_times == validation_start) + 1L])

  # Each fit uses only outcomes that would have been known when it was made.
  train <- data[data$opened < validation_start & data$closed < validation_start, ]
  validation <- data[data$opened >= validation_start & data$opened < test_start &
                       data$closed < test_start, ]
  final_train <- data[data$opened < test_start & data$closed < test_start, ]
  test <- data[data$opened >= test_start, ]
  if (min(nrow(train), nrow(validation), nrow(test)) < 2L) {
    stop("Too few cases in a split. Use a sample covering a longer time span.")
  }

  # Learn the categories from the fitting rows. One indicator per observed
  # category lets ridge shrink all categories symmetrically.
  make_design <- function(fit, assess, predictors) {
    fit_columns <- list()
    assess_columns <- list()
    for (variable in predictors) {
      clean <- function(x) {
        x <- trimws(as.character(x))
        x[is.na(x) | x == ""] <- "(Missing)"
        x
      }

      x <- clean(fit[[variable]])
      new_x <- clean(assess[[variable]])
      counts <- table(x)
      keep <- names(counts)
      if (variable %in% c("service_name", "analysis_neighborhood")) {
        keep <- union(names(counts[counts >= rare_cutoff]), "(Missing)")
        x[!x %in% keep] <- "(Other)"
        new_x[!new_x %in% keep] <- "(Other)"
      }

      # A category absent from the fitting rows has no fitted effect.
      levels_seen <- sort(unique(x))
      fit_columns[[variable]] <- outer(x, levels_seen, "==") * 1
      assess_columns[[variable]] <- outer(new_x, levels_seen, "==") * 1
    }
    list(fit = do.call(cbind, fit_columns),
         assess = do.call(cbind, assess_columns))
  }

  predict_ridge <- function(design, y, penalty) {
    # Center the columns to fit an intercept without penalizing it.
    # Minimize mean squared fitting error plus penalty * sum(beta^2).
    column_means <- colMeans(design$fit)
    x <- sweep(design$fit, 2L, column_means)
    new_x <- sweep(design$assess, 2L, column_means)
    y_mean <- mean(y)
    coefficients <- solve(crossprod(x) + diag(nrow(x) * penalty, ncol(x)),
                          crossprod(x, y - y_mean))
    pmax(0, drop(y_mean + new_x %*% coefficients))
  }

  penalties <- 10^seq(-5, 2, by = 0.5)
  options <- c("No location" = "", "Police district" = location[1L],
               "Supervisor district" = location[2L], "Neighborhood" = location[3L])
  rows <- lapply(names(options), function(model) {
    predictors <- c(baseline, unname(options[model]))
    predictors <- predictors[nzchar(predictors)]
    validation_design <- make_design(train, validation, predictors)
    validation_mae <- vapply(penalties, function(penalty) {
      predicted <- predict_ridge(validation_design, train$hours_to_close, penalty)
      mean(abs(validation$hours_to_close - predicted))
    }, numeric(1))
    best <- which.min(validation_mae)

    test_design <- make_design(final_train, test, predictors)
    predicted <- predict_ridge(test_design, final_train$hours_to_close, penalties[best])
    error <- test$hours_to_close - predicted
    data.frame(model = model, penalty = penalties[best],
               validation_mae_hours = validation_mae[best],
               test_mae_hours = mean(abs(error)),
               test_rmse_hours = sqrt(mean(error^2)))
  })

  results <- do.call(rbind, rows)
  baseline_mae <- results$test_mae_hours[results$model == "No location"]
  results$mae_improvement_minutes <- 60 * (baseline_mae - results$test_mae_hours)

  cat("Validation starts:", format(validation_start), "\n")
  cat("Test starts:", format(test_start), "\n")
  cat("Initial train / validation / refit train / test rows:",
      nrow(train), "/", nrow(validation), "/", nrow(final_train), "/", nrow(test), "\n")
  results[order(results$test_mae_hours), ]
}

location_results <- compare_location(cases_clean)
print(location_results, digits = 4, row.names = FALSE)
