### Helper functions to run the scripts ----

## trauma_rate.R script helper ----
# fit & bootstrap per-subset function 
fit_region_spline <- function(df_region, pred_x = NULL, n_boot = 500, spar = 0.85, min_unique_x = 4) {
  # df_region: data.frame with Mid.date and prop_smooth
  # returns: list with data.frame 'pred_df' (x, fit_pct, ci_low, ci_high) and counts/n_used
  
  df_region <- df_region |> 
    arrange(Mid.date)
  
  n_rows <- nrow(df_region)
  n_unique_x <- length(unique(df_region$Mid.date))
  
  if (n_rows < 3 || n_unique_x < min_unique_x) {
    warning("Insufficient data for spline (rows=", n_rows, ", unique_x=", n_unique_x, "). Skipping fit.")
    return(NULL)
  }
  
  # logit transform
  df_region <- df_region |> 
    mutate(logit_prop = qlogis(prop_smooth))
  
  if (any(!is.finite(df_region$logit_prop))) {
    warning("Non-finite logit_prop encountered; skipping region.")
    return(NULL)
  }
  
  # Primary fit
  fit_primary <- tryCatch(
    smooth.spline(x = df_region$Mid.date, y = df_region$logit_prop, spar = spar),
    error = function(e) { warning("smooth.spline error: ", e$message); NULL }
  )
  if (is.null(fit_primary)) return(NULL)
  
  # Use specified pred_x or fit_primary$x
  pred_x_use <- if (is.null(pred_x)) fit_primary$x else pred_x
  pred_y_logit <- predict(fit_primary, x = pred_x_use)$y
  pred_y_pct <- plogis(pred_y_logit) * 100
  
  # Bootstrap
  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(pred_x_use))
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(nrow(df_region)), replace = TRUE)
    samp <- df_region[idx, ]
    if (length(unique(samp$Mid.date)) < min_unique_x) {
      boot_mat[b, ] <- NA_real_
      next
    }
    bfit <- tryCatch(smooth.spline(x = samp$Mid.date, y = samp$logit_prop, spar = spar),
                     error = function(e) NULL, warning = function(w) NULL)
    if (is.null(bfit)) {
      boot_mat[b, ] <- NA_real_
      next
    }
    preds_logit <- tryCatch(predict(bfit, x = pred_x_use)$y, error = function(e) rep(NA_real_, length(pred_x_use)))
    boot_mat[b, ] <- plogis(preds_logit) * 100
  }
  valid_boots <- boot_mat[rowSums(is.na(boot_mat)) < ncol(boot_mat), , drop = FALSE]
  if (nrow(valid_boots) < max(10, round(n_boot * 0.05))) {
    warning("Few valid bootstrap replicates for region (", nrow(valid_boots), "). CIs may be unreliable.")
  }
  ci_low <- apply(valid_boots, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_high <- apply(valid_boots, 2, quantile, probs = 0.975, na.rm = TRUE)
  # replace any NA CI positions with fit value
  nas <- which(!is.finite(ci_low) | !is.finite(ci_high))
  if (length(nas) > 0) {
    ci_low[nas] <- pred_y_pct[nas]
    ci_high[nas] <- pred_y_pct[nas]
  }
  
  pred_df <- data.frame(x = pred_x_use, fit_pct = pred_y_pct, ci_low = ci_low, ci_high = ci_high)
  attr(pred_df, "n_rows") <- n_rows
  attr(pred_df, "n_valid_boot") <- nrow(valid_boots)
  return(pred_df)
}

## ethnographic_comparison.R helper ----
# Vargha-Delaney A
vda_A <- function(x, y) {
  r <- rank(c(x, y))
  R1 <- sum(r[seq_along(x)])
  (R1 / length(x) - (length(x) + 1) / 2) / length(y)
}

# bootstrap median difference between arrows
boot_median_diff <- function(x, y, n_boot) {
  replicate(n_boot, median(sample(x, length(x), replace = TRUE)) -
              median(sample(y, length(y), replace = TRUE)))
}

## classifying.R helper ----
# apply posterior_epred and get the median
post_median <- function(fit, newdata) {
  apply(posterior_epred(fit, newdata = newdata), 2, median)
}
