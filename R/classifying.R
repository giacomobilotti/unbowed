### Classifying Palpa points ----
# This script was originally created by E. Marsh 

# Classify Palpa projectile points as dart or arrow against hafted archaeological points
# Script written by Claude AI (model 4.7). 
# It was checked and processed in R Studio by E. Marsh.
# Prior to publicaiton, G. Bilotti re run the script and integrated it into the repository

# Bayesian logistic regression on log-transformed metric attributes, trained
# on the Appendix A hafted-point reference compilation (Marsh et al. 2024).
# Classifies all Palpa projectile points with valid width and thickness
# measurements (complete and incomplete), under five model specifications:
# one bivariate (width + thickness) and four univariate (width, thickness,
# TCSA, TCSP). Each Palpa point receives a posterior median P(dart) for each
# model.
#
#
# Incomplete points (Section 8 of this script):
#   Fragment measurements are lower bounds on the original dimensions, since
#   breakage removes material. Each fragment is classified using only the
#   dimension(s) preserved through breakage. We diagnose which dimension was
#   lost by comparing log(width / thickness) of the fragment to the empirical
#   2.5th–97.5th percentile range of log(width / thickness) among complete
#   Palpa points (the Palpa-specific allometric reference):
#     - ratio below the range → width disproportionately reduced; classify
#                                with the thickness-only model.
#     - ratio above the range → thickness disproportionately reduced;
#                                classify with the width-only model.
#     - ratio within the range → length lost (tip/base); both width and
#                                 thickness intact; classify with the
#                                 bivariate model.
#   The resulting assigned P(dart) is stored as P_dart_assigned. A fragment
#   classified as DART under this scheme is robust (the true original was
#   at least as large and would classify the same way or more strongly
#   dartward). A fragment classified as ARROW is less certain.
#
# Rounding-error propagation (Section 9):
#   Palpa measurements were recorded to the nearest 1 mm, so the
#   true value of each recorded width or thickness is uniformly distributed
#   in [recorded − 0.5, recorded + 0.5] mm.  We propagate the Palpa rounding uncertainty
#   through to P(dart) by Monte
#   Carlo: M = 100 perturbed copies per complete point, with width and
#   thickness independently shifted by Uniform(−0.5, 0.5); TCSA and TCSP
#   are recomputed from the perturbed values.
#
# Inputs (current working directory):
#   - Hafted_points_Appendix_A.xlsx   (Marsh et al. 2024)
#   - Palpa_points.xlsx               (this study)
#
# Outputs:
#   - Palpa_classified.xlsx           (P(dart) under five models, all points)
#   - brms_fits.rds                   (fitted models)
#   - Palpa_classification_plot.svg   (bivariate decision surface)
#   - Palpa_pdart_summary.svg         (P(dart) across the five models)
#   - Palpa_width_logistic_plot.svg   (univariate width, Haas and Kelly-comparable)
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Setup
# -----------------------------------------------------------------------------

library(readxl)
library(writexl)
library(brms)
library(dplyr)
library(ggplot2)
library(ggnewscale)  # multiple fill / shape scales in the same plot

# Fix the random-number generator's starting state. R's default RNG produces
# a deterministic sequence of pseudorandom numbers once the seed is set, so
# every call to runif() (rounding-error perturbation, Section 9) and to
# brms/Stan (MCMC chains, Section 5) returns the same draws on every run.
# Without this line the results would be very close from run to run — but
# not bit-identical — and supplementary outputs (perturbed predictions,
# posterior medians at the 4th decimal place) would shift slightly. The
# specific integer (today's date in YYYYMMDD form) is arbitrary; any fixed
# value would do.

set.seed(20260518)

# Palpa integrity-grade palette

palpa_fill <- c("A+B" = "#a7ea52",  # lime
                "C"   = "#34ac8b",  # teal
                "D"   = "#beebdf")  # mint

# Diverging gradient for P(dart): red (arrow) → cream (boundary) → blue (dart).

pdart_low  <- "#d7191c"  # P(dart) = 0   → arrow
pdart_mid  <- "#fff5cc"  # P(dart) = 0.5 → boundary
pdart_high <- "#2c7bb6"  # P(dart) = 1   → dart

# Palpa point shapes: filled square = complete, filled circle = incomplete.

palpa_shape <- c("Complete" = 22, "Incomplete" = 21)


# -----------------------------------------------------------------------------
# 2. Load and prepare training data
#
# Source: Appendix A of Marsh et al. (2024), a hafted-point compilation
# spanning the Americas. We keep only specimens unambiguously identified as
# Dart or Arrow (excluding “Unclear” cases). 
# Training is restricted to complete cases on the four shared predictors
# (thickness, width, TCSA, TCSP), since the bivariate and TCSA/TCSP univariate
# models require these to be present.
# -----------------------------------------------------------------------------

haft <- read_excel("Hafted_points_Appendix_A.xlsx", sheet = "Data table")

train <- haft %>%
  filter(`Weapon Type` %in% c("Dart", "Arrow")) %>%
  transmute(
    is_dart   = as.integer(`Weapon Type` == "Dart"),
    thickness = `Point thickness (mm)`,
    width     = `Maximum point width (mm)`,
    tcsa      = `Tip cross-section area (mm²)`,
    tcsp      = `Tip cross-section perimeter (mm)`
  ) %>%
  filter(complete.cases(.))

cat("\n--- Training set composition ---\n")
cat("Training n      :", nrow(train), "\n")
cat("  Darts         :", sum(train$is_dart), "\n")
cat("  Arrows        :", sum(!train$is_dart), "\n")
cat("  Ratio dart:arrow :",
    sprintf("%.2f : 1", sum(train$is_dart) / sum(!train$is_dart)), "\n")
cat("Compare Buchanan et al. (2025): 51 darts, 220 arrows (1 : 4.31)\n")


# -----------------------------------------------------------------------------
# 3. Diagnostics: predictor correlation
#
# Width and thickness are physically linked (a wider point is on average
# thicker). Haas and Kelly (2026) report a strong correlation between width
# and thickness in the Thomas (1978) / Shott (1997) reference (t ≈ 11.9,
# p < 0.01) and recommend treating width and thickness as collinear. We
# report Pearson r and the implied VIF for the bivariate model on log-
# transformed predictors.
# -----------------------------------------------------------------------------

r_wt   <- cor(log(train$width), log(train$thickness))
vif_wt <- 1 / (1 - r_wt^2)

cat("\n--- Predictor correlation ---\n")
cat("Pearson r [log(width), log(thickness)] :", round(r_wt, 3), "\n")
cat("VIF for bivariate model               :", round(vif_wt, 2), "\n")
cat("Rule of thumb: VIF > 5 indicates appreciable multicollinearity.\n")


# -----------------------------------------------------------------------------
# 4. Load and prepare Palpa points
#
# Includes both Complete and Incomplete points. Incomplete-point measurements
# are lower bounds on the original dimensions (breakage removes material);
# they are still classified by the same models and flagged in the output for
# interpretation.
#
# TCSA and TCSP in Palpa_points.xlsx are computed as ½·w·t and
# 2·√(w² + t²) respectively (mm² and mm), matching the Appendix A
# convention, so they are loaded as-is.
# -----------------------------------------------------------------------------

palpa <- read_excel("Palpa_points.xlsx", sheet = "1. Projectile points")

palpa_classified <- palpa %>%
  rename(
    thickness    = `Max thickness (mm)`,
    width        = `Max width (mm)`,
    tcsa         = `Tip cross-section area (mm²)`,
    tcsp         = `Tip cross-section perimeter (mm)`,
    completeness = `Complete? (measurements of fragments are shaded since they may be unreliable for metric comparisons)`
  ) %>%
  filter(completeness %in% c("Complete", "Incomplete"),
         !is.na(thickness), !is.na(width)) %>%
  mutate(
    completeness = factor(completeness, levels = c("Complete", "Incomplete")),
    is_complete  = completeness == "Complete"
  )

cat("\n--- Palpa points ---\n")
cat("Total classifiable Palpa points n =", nrow(palpa_classified), "\n")
cat("  Complete   :", sum(palpa_classified$is_complete), "\n")
cat("  Incomplete :", sum(!palpa_classified$is_complete), "\n")


# -----------------------------------------------------------------------------
# 5. Fit five Bayesian logistic regression models
#
# Weakly informative priors: N(0, 5) on the intercept and N(0, 2.5) on slopes.
# On the log-odds scale, N(0, 5) puts prior mass essentially uniformly across
# P(dart) ∈ (0.01, 0.99) centred on 0.5 — equivalent to the vague 50:50 prior
# recommended by Haas and Kelly (2026). The slope prior is weakly informative
# around zero.
#
# brms uses Stan via the No-U-Turn Sampler. Defaults: 4 chains, 2,000
# iterations each (1,000 warmup, 1,000 post-warmup), total 4,000 posterior
# draws.
# -----------------------------------------------------------------------------

priors <- c(
  prior(normal(0, 5),   class = "Intercept"),
  prior(normal(0, 2.5), class = "b")
)

cat("\n--- Fitting brms models (5 fits, ~3–5 min total) ---\n")

fit_bivar <- brm(is_dart ~ log(thickness) + log(width),
                 data = train, family = bernoulli(),
                 prior = priors, seed = 101, refresh = 0, silent = 2)

fit_thick <- brm(is_dart ~ log(thickness),
                 data = train, family = bernoulli(),
                 prior = priors, seed = 102, refresh = 0, silent = 2)

fit_width <- brm(is_dart ~ log(width),
                 data = train, family = bernoulli(),
                 prior = priors, seed = 103, refresh = 0, silent = 2)

fit_tcsa  <- brm(is_dart ~ log(tcsa),
                 data = train, family = bernoulli(),
                 prior = priors, seed = 104, refresh = 0, silent = 2)

fit_tcsp  <- brm(is_dart ~ log(tcsp),
                 data = train, family = bernoulli(),
                 prior = priors, seed = 105, refresh = 0, silent = 2)

# In-sample classification accuracy on the training set, by model.
# Strictly a sanity check; not a substitute for proper cross-validation.
acc <- function(fit) {
  p <- apply(posterior_epred(fit), 2, median)
  mean(round(p) == train$is_dart)
}
cat("\n--- Training-set in-sample accuracy ---\n")
cat("  bivariate :", round(acc(fit_bivar), 3), "\n")
cat("  thickness :", round(acc(fit_thick), 3), "\n")
cat("  width     :", round(acc(fit_width), 3), "\n")
cat("  TCSA      :", round(acc(fit_tcsa),  3), "\n")
cat("  TCSP      :", round(acc(fit_tcsp),  3), "\n")


# -----------------------------------------------------------------------------
# 6. Predict posterior median P(dart) for every Palpa point
#
# Applies to both complete and incomplete points. Per-point credible
# intervals are available from the same posterior_epred() output if needed.
# -----------------------------------------------------------------------------

post_median <- function(fit, newdata) {
  apply(posterior_epred(fit, newdata = newdata), 2, median)
}

palpa_classified$P_dart_bivariate <- post_median(fit_bivar, palpa_classified)
palpa_classified$P_dart_thickness <- post_median(fit_thick, palpa_classified)
palpa_classified$P_dart_width     <- post_median(fit_width, palpa_classified)
palpa_classified$P_dart_TCSA      <- post_median(fit_tcsa,  palpa_classified)
palpa_classified$P_dart_TCSP      <- post_median(fit_tcsp,  palpa_classified)


# -----------------------------------------------------------------------------
# 7. Width threshold μ (univariate width model)
#
# Method: for every posterior draw, find the width at which the linear
# predictor crosses zero (i.e., where P(dart) = 0.5). This yields a posterior
# distribution of the threshold from the brms univariate width fit.
#
# -----------------------------------------------------------------------------

grid_w <- data.frame(width = seq(5, 40, length.out = 2000))
lp     <- posterior_linpred(fit_width, newdata = grid_w)  # draws × widths

threshold_draws <- apply(lp, 1, function(eta) {
  approx(eta, grid_w$width, xout = 0)$y
})
threshold_draws <- threshold_draws[is.finite(threshold_draws)]

threshold_med <- median(threshold_draws)              # non-parametric point estimate
threshold_ci  <- quantile(threshold_draws,            # non-parametric interval
                          c(0.025, 0.975))

cat("\n--- Width threshold (P(dart) = 0.5) ---\n")
cat(sprintf("  This study (Marsh et al. 2024 reference): %.2f mm",
            threshold_med),
    sprintf(" (95%% CrI [posterior percentiles]: %.2f – %.2f mm)\n",
            threshold_ci[1], threshold_ci[2]))
cat("  Haas and Kelly (2026, Buchanan et al. (2025) expanded reference): ",
    "16.6 – 18.5 mm (median 17.6 mm)\n", sep = "")


# -----------------------------------------------------------------------------
# 8. Breakage diagnosis and assigned-model P(dart) for incomplete points
#
# Fragment measurements are lower bounds on the original dimensions. To
# classify incomplete points honestly, we diagnose which dimension was lost
# in breakage, then classify using only the trustworthy dimension(s).
#
# Reference distribution: log(width / thickness) for the complete Palpa
# points. The empirical 2.5th–97.5th percentile interval of that ratio
# defines the "normal" range. An incomplete point whose ratio falls:
#   - below the lower bound  → width was disproportionately reduced
#                              (side broken; thickness intact); use the
#                              thickness-only model.
#   - above the upper bound  → thickness was disproportionately reduced
#                              (face spalled; width intact); use the
#                              width-only model. This pattern is rare in
#                              lithic reduction but included for
#                              completeness.
#   - within the bounds      → ratio is plausible; breakage is most likely
#                              along the length (tip or base missing) and
#                              both width and thickness are intact; use
#                              the bivariate model.
#
# This diagnostic uses the complete Palpa points as the reference (not the
# Appendix A hafted set), because the *Palpa-specific* width/thickness
# allometry is what the local breakage diagnosis needs to compare against.
# -----------------------------------------------------------------------------

complete_subset <- palpa_classified %>% filter(is_complete)
complete_subset$log_ratio <- log(complete_subset$width /
                                   complete_subset$thickness)

ref_lo  <- as.numeric(quantile(complete_subset$log_ratio, 0.025))
ref_hi  <- as.numeric(quantile(complete_subset$log_ratio, 0.975))
ref_med <- median(complete_subset$log_ratio)

cat("\n--- Breakage-diagnosis reference (complete Palpa points) ---\n")
cat(sprintf("  log(width/thickness): median %.3f, 95%% range [%.3f, %.3f]\n",
            ref_med, ref_lo, ref_hi))

# Diagnose every Palpa point (complete points are diagnosed as "intact",
# incomplete points get a breakage_pattern label).
palpa_classified <- palpa_classified %>%
  mutate(
    log_ratio        = log(width / thickness),
    breakage_pattern = case_when(
      !is_complete & log_ratio < ref_lo ~ "width_broken",
      !is_complete & log_ratio > ref_hi ~ "thickness_broken",
      !is_complete                      ~ "length_broken",
      TRUE                              ~ "intact"
    )
  )

cat("\n--- Breakage diagnosis for incomplete points ---\n")
print(table(palpa_classified$breakage_pattern[!palpa_classified$is_complete]))

# Assigned-model P(dart): pick the model that uses only the trustworthy
# dimension(s). Complete points get the bivariate model by default.
palpa_classified$P_dart_assigned <- with(palpa_classified, dplyr::case_when(
  breakage_pattern == "width_broken"     ~ P_dart_thickness,
  breakage_pattern == "thickness_broken" ~ P_dart_width,
  TRUE                                   ~ P_dart_bivariate  # intact or length_broken
))
palpa_classified$model_used <- with(palpa_classified, dplyr::case_when(
  breakage_pattern == "width_broken"     ~ "thickness_only",
  breakage_pattern == "thickness_broken" ~ "width_only",
  TRUE                                   ~ "bivariate"
))


# -----------------------------------------------------------------------------
# 9. Measurement-error propagation: rounding-induced uncertainty in P(dart)
#
# Palpa measurements were recorded to the nearest 1 mm.
#
# We propagate this rounding uncertainty through to P(dart) by Monte Carlo.
# For each complete Palpa point, M = 100 perturbed copies are generated:
# width and thickness are each shifted independently by Uniform(−0.5, 0.5).
# TCSA = ½·w·t and TCSP = 2·√(w² + t²) are then RECOMPUTED from the
# perturbed (w, t) values (they are derived measurements and inherit
# uncertainty from width and thickness; perturbing them independently would
# be incorrect). Each of the five brms models then predicts P(dart) for
# every perturbed copy. The 2.5th and 97.5th percentiles per (point ×
# model) give a non-parametric 95% interval on P(dart) attributable to
# rounding alone.
#
# Note on the relative widths of error bars: TCSA and TCSP intervals are
# typically wider than width-only or thickness-only intervals. This is
# expected, not a bug: a ±0.5 mm shift on a 5 mm thickness is a 10% change
# in thickness; TCSA inherits proportional uncertainty from both factors,
# so a small point's TCSA can move by 15–20% under perturbation. The
# logistic regression then maps that into a wider P(dart) range, especially
# for points sitting near the decision boundary. This propagation finding
# is itself informative — it shows that TCSA/TCSP-based classifications are
# more sensitive to rounded measurements than width- or thickness-based
# ones, particularly for small points.
#
# Incomplete points are not included here because their measurement
# uncertainty is dominated by breakage (handled in Section 8), not calliper
# precision.
# -----------------------------------------------------------------------------

M <- 100  # Monte Carlo perturbations per point

complete_idx  <- which(palpa_classified$is_complete)
n_complete    <- length(complete_idx)

expanded <- palpa_classified[rep(complete_idx, each = M), , drop = FALSE]
expanded$row_id <- rep(complete_idx, each = M)

expanded$width     <- expanded$width     + runif(n_complete * M, -0.5, 0.5)
expanded$thickness <- expanded$thickness + runif(n_complete * M, -0.5, 0.5)
expanded$tcsa      <- 0.5 * expanded$width * expanded$thickness
expanded$tcsp      <- 2   * sqrt(expanded$width^2 + expanded$thickness^2)

cat(sprintf("\n--- Rounding-error Monte Carlo (M = %d, n = %d complete points) ---\n",
            M, n_complete))
cat("  Running posterior_epred() under perturbed measurements...\n")

# Predict per model from the existing brms fits. apply() collapses the
# (draws × rows) posterior_epred() output to a per-row posterior median;
# we then aggregate those medians within each Palpa point's M perturbations
# to get the rounding-induced 95% interval.
rounding_pred <- function(fit, newdata) {
  apply(posterior_epred(fit, newdata = newdata), 2, median)
}

expanded$p_bivar <- rounding_pred(fit_bivar, expanded)
expanded$p_thick <- rounding_pred(fit_thick, expanded)
expanded$p_width <- rounding_pred(fit_width, expanded)
expanded$p_tcsa  <- rounding_pred(fit_tcsa,  expanded)
expanded$p_tcsp  <- rounding_pred(fit_tcsp,  expanded)

# Per Palpa point, summarise the M perturbations: lo / median / hi for
# each model.
rounding_summary <- expanded %>%
  group_by(row_id) %>%
  summarise(
    P_dart_bivariate_lo = quantile(p_bivar, 0.025),
    P_dart_bivariate_hi = quantile(p_bivar, 0.975),
    P_dart_thickness_lo = quantile(p_thick, 0.025),
    P_dart_thickness_hi = quantile(p_thick, 0.975),
    P_dart_width_lo     = quantile(p_width, 0.025),
    P_dart_width_hi     = quantile(p_width, 0.975),
    P_dart_TCSA_lo      = quantile(p_tcsa,  0.025),
    P_dart_TCSA_hi      = quantile(p_tcsa,  0.975),
    P_dart_TCSP_lo      = quantile(p_tcsp,  0.025),
    P_dart_TCSP_hi      = quantile(p_tcsp,  0.975),
    .groups = "drop"
  )

# Merge rounding-error bounds back onto the main table. Incomplete points stay NA
# in the *_lo / *_hi columns. `dplyr::` qualifier prevents collision with
# MASS::select (which loads as a brms/Stan dependency on some systems).
palpa_classified <- palpa_classified %>%
  mutate(row_id = dplyr::row_number()) %>%
  dplyr::left_join(rounding_summary, by = "row_id") %>%
  dplyr::select(-row_id)

cat("  Done. Rounding-error 95% intervals added to output.\n")


# -----------------------------------------------------------------------------
# 10. Save outputs
# -----------------------------------------------------------------------------

write_xlsx(palpa_classified, "Palpa_classified.xlsx")

saveRDS(list(bivar = fit_bivar, thick = fit_thick, width = fit_width,
             tcsa  = fit_tcsa,  tcsp  = fit_tcsp, train = train),
        "brms_fits.rds")

cat("\n--- Saved ---\n")
cat("  Palpa_classified.xlsx\n")
cat("  brms_fits.rds\n")


# -----------------------------------------------------------------------------
# 11. Plotting setup shared across all three plots
# -----------------------------------------------------------------------------

# Add integrity-grade grouping (A+B / C / D), used as fill color in every
# plot. Points with a missing integrity grade are dropped from plotting.
# `plot_palpa` is the full set including fragments; `plot_palpa_complete`
# is restricted to complete points only and is used by Plots 1 and 3.
plot_palpa <- palpa_classified %>%
  mutate(grade_group = case_when(
    `Integrity grade` %in% c("A", "B") ~ "A+B",
    `Integrity grade` == "C"           ~ "C",
    `Integrity grade` == "D"           ~ "D",
    TRUE                               ~ NA_character_
  )) %>%
  filter(!is.na(grade_group))

plot_palpa_complete <- plot_palpa %>% filter(is_complete)

train_lab <- train %>%
  mutate(weapon = ifelse(is_dart == 1L, "Dart", "Arrow"))


# -----------------------------------------------------------------------------
# 12. Plot 1: bivariate decision surface (width × thickness; Figure 7)
#
# Background: posterior median P(dart) from the bivariate model on a grid of
# (width, thickness). Foreground: hafted reference points (open shapes) and
# complete Palpa points as filled squares colored by integrity grade.
# Incomplete points are excluded because their width and thickness
# measurements are lower bounds on the original dimensions and would plot
# in misleading positions.
# Width on x-axis.
# -----------------------------------------------------------------------------

w_range <- range(c(plot_palpa_complete$width, train_lab$width),         na.rm = TRUE)
t_range <- range(c(plot_palpa_complete$thickness, train_lab$thickness), na.rm = TRUE)

grid <- expand.grid(
  width     = seq(w_range[1] - 1,   w_range[2] + 2,   length.out = 150),
  thickness = seq(t_range[1] - 0.3, t_range[2] + 0.5, length.out = 150)
)
grid$p_dart <- apply(posterior_epred(fit_bivar, newdata = grid), 2, median)

p1 <- ggplot() +
  # Posterior median P(dart) surface
  geom_raster(data = grid,
              aes(width, thickness, fill = p_dart),
              alpha = 0.55, interpolate = TRUE) +
  scale_fill_gradient2(
    low = pdart_low, mid = pdart_mid, high = pdart_high,
    midpoint = 0.5, limits = c(0, 1),
    name = "P(dart)\nbivariate"
  ) +
  # Decision boundary
  geom_contour(data = grid,
               aes(width, thickness, z = p_dart),
               breaks = 0.5, colour = "black", linewidth = 0.4) +
  # Hafted reference, faded open shapes
  geom_point(data = train_lab,
             aes(width, thickness, shape = weapon),
             colour = "gray25", fill = NA,
             alpha = 0.45, size = 2.2, stroke = 0.7) +
  scale_shape_manual(values = c("Arrow" = 1, "Dart" = 2),
                     name = "Hafted reference") +
  # Separate fill scale for Palpa points (complete only; all squares).
  new_scale_fill() +
  geom_point(data = plot_palpa_complete,
             aes(width, thickness, fill = grade_group),
             shape = 22, colour = "black", stroke = 0.5,
             size = 2.8, alpha = 0.95) +
  scale_fill_manual(
    values = palpa_fill, name = "Palpa integrity",
    guide = guide_legend(
      override.aes = list(shape = 22, colour = "black", size = 3.2)
    )
  ) +
  labs(
    x = "Maximum width (mm)",
    y = "Maximum thickness (mm)",
    title = "Palpa points against hafted dart/arrow reference",
    subtitle = sprintf(
      "Background: posterior median P(dart) from bivariate Bayesian logistic regression. n = %d complete points.",
      nrow(plot_palpa_complete))
  ) +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.box       = "vertical")

ggsave("Palpa_classification_plot.svg", p1,
       width = 10, height = 7, device = "svg")


# -----------------------------------------------------------------------------
# 13. Plot 2: P(dart) summary with rounding error and incomplete points (Figure S9)
#
# For complete Palpa points, plot the posterior median P(dart) under each of
# the five models in five separate rows, with horizontal error bars showing
# the 95% interval from ±0.5 mm rounding (Section 9). A sixth row shows
# incomplete points using the assigned-model P(dart) from Section 8 — i.e.,
# each fragment is classified using only the dimension(s) preserved through
# breakage (bivariate, thickness-only, or width-only depending on the
# breakage_pattern). Rounding error is not shown for incomplete points
# because their measurement uncertainty is dominated by breakage, not
# calliper precision.
#
# Lets the reader assess:
#   (i)   classification robustness across model choices,
#   (ii)  classification confidence (distance from P = 0.5),
#   (iii) sensitivity to ±0.5 mm rounding error (length of error bars),
#   (iv)  integrity-grade pattern (do A+B, C, D points behave alike?),
#   (v)   how fragmented points compare to complete ones.
# -----------------------------------------------------------------------------

models <- c("bivariate", "thickness", "width", "TCSA", "TCSP")

# Subsets used to label rows. complete_plot drives rows 1–5; incomplete_plot
# drives row 6 (assigned model).
complete_plot   <- plot_palpa %>% filter(is_complete)
incomplete_plot <- plot_palpa %>% filter(!is_complete)

n_complete   <- nrow(complete_plot)
n_incomplete <- nrow(incomplete_plot)

# Row labels for the y-axis. Short keys (`model_key`) preserve the column
# mapping used below to pull the right `P_dart_*` columns; long-form
# `row_label` is shown on the plot.
row_labels <- c(
  bivariate  = sprintf("bivariate\n(width + thickness, n = %d)", n_complete),
  thickness  = sprintf("thickness only\n(n = %d)", n_complete),
  width      = sprintf("width only\n(n = %d)", n_complete),
  TCSA       = sprintf("TCSA\n(n = %d)", n_complete),
  TCSP       = sprintf("TCSP\n(n = %d)", n_complete),
  assigned   = sprintf("incomplete points\n(model varies by breakage, n = %d)",
                       n_incomplete)
)
row_order <- c("bivariate", "thickness", "width", "TCSA", "TCSP", "assigned")

palpa_long_complete <- do.call(rbind, lapply(models, function(m) {
  data.frame(
    grade_group  = complete_plot$grade_group,
    completeness = factor("Complete", levels = c("Complete", "Incomplete")),
    model_key    = factor(m, levels = row_order),
    p_dart       = complete_plot[[paste0("P_dart_", m)]],
    p_lo         = complete_plot[[paste0("P_dart_", m, "_lo")]],
    p_hi         = complete_plot[[paste0("P_dart_", m, "_hi")]]
  )
}))

# Sixth row: incomplete points using the assigned model. No rounding bars.
palpa_long_incomplete <- data.frame(
  grade_group  = incomplete_plot$grade_group,
  completeness = factor("Incomplete", levels = c("Complete", "Incomplete")),
  model_key    = factor("assigned", levels = row_order),
  p_dart       = incomplete_plot$P_dart_assigned,
  p_lo         = NA_real_,
  p_hi         = NA_real_
)

palpa_long <- rbind(palpa_long_complete, palpa_long_incomplete)

# Per-row dart-count summary text.
summary_df <- palpa_long %>%
  group_by(model_key) %>%
  summarise(n      = n(),
            n_dart = sum(p_dart >= 0.5),
            pct    = round(100 * n_dart / n),
            .groups = "drop") %>%
  mutate(label = sprintf("P(dart) ≥ 0.5: %d / %d (%d%%)",
                         n_dart, n, pct))

# Background P(dart) gradient: 200 narrow vertical bands across x ∈ [0, 1].
n_strip <- 200
xs <- seq(0, 1, length.out = n_strip + 1)
bg_strip <- data.frame(
  xmin = xs[-(n_strip + 1)],
  xmax = xs[-1],
  x    = (xs[-(n_strip + 1)] + xs[-1]) / 2,
  ymin = 0.5,
  ymax = length(row_order) + 0.5
)

# Per-row vertical jitter offset, deterministic so error bars line up with
# their points.
set.seed(7)
palpa_long$y_offset <- runif(nrow(palpa_long), -0.22, 0.22)
palpa_long$y_pos    <- as.numeric(palpa_long$model_key) + palpa_long$y_offset

p2 <- ggplot() +
  # Background P(dart) gradient
  geom_rect(data = bg_strip,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill = x),
            alpha = 0.22, inherit.aes = FALSE) +
  scale_fill_gradient2(
    low = pdart_low, mid = pdart_mid, high = pdart_high,
    midpoint = 0.5, limits = c(0, 1),
    name = NULL, guide = "none"
  ) +
  # Decision boundary
  geom_vline(xintercept = 0.5, linetype = "dashed",
             colour = "gray25", linewidth = 0.5) +
  # Horizontal separator between five-model block and incomplete row
  geom_hline(yintercept = length(models) + 0.5,
             colour = "gray60", linetype = "dotted", linewidth = 0.4) +
  # Rounding-error 95% intervals (NA-omitted automatically for incomplete row)
  geom_segment(data = palpa_long %>% filter(!is.na(p_lo)),
               aes(x = p_lo, xend = p_hi, y = y_pos, yend = y_pos),
               colour = "gray35", alpha = 0.45, linewidth = 0.35,
               inherit.aes = FALSE) +
  # Separate fill and shape scales for Palpa points
  new_scale_fill() +
  geom_point(data = palpa_long,
             aes(x = p_dart, y = y_pos,
                 fill = grade_group, shape = completeness),
             colour = "black", stroke = 0.4,
             size = 2.2, alpha = 0.88) +
  scale_fill_manual(
    values = palpa_fill, name = "Palpa integrity",
    # Force the legend keys to render as filled squares with black outline,
    # otherwise ggplot draws them as the default circle (no fill aesthetic),
    # which renders as a black dot.
    guide = guide_legend(
      override.aes = list(shape = 22, colour = "black", size = 3.2)
    )
  ) +
  scale_shape_manual(
    values = palpa_shape, name = "Palpa point",
    # The Palpa-point legend shows the actual shapes (square/circle) with
    # a neutral fill so the *shape* contrast is what reads, not the colour.
    guide = guide_legend(
      override.aes = list(fill = "gray70", colour = "black", size = 3.2)
    )
  ) +
  # Per-row dart-count text, above each row.
  geom_text(data = summary_df,
            aes(x = 0.02,
                y = as.numeric(model_key),
                label = label),
            hjust = 0, vjust = -1.7, size = 3.2, colour = "gray15",
            inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = 0)) +
  scale_y_reverse(
    breaks = seq_along(row_order),
    labels = row_labels[row_order],
    expand = expansion(add = 0.5)
  ) +
  labs(
    x = "Posterior median P(dart)",
    y = NULL,
    title    = "P(dart) summary across five models for Palpa points",
    subtitle = sprintf(
      paste0("%d complete points (rows 1–5; error bars = 95%% interval from ±0.5 mm measurement rounding) ",
             "+ %d incomplete points (row 6; assigned model varies by breakage diagnosis)."),
      n_complete, n_incomplete)
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(lineheight = 0.9))

ggsave("Palpa_pdart_summary.svg", p2,
       width = 11, height = 7, device = "svg")


# -----------------------------------------------------------------------------
# 14. Plot 3: univariate width logistic curve (Figure 8)
#
#Follows Haas and Kelly (2026, figure 3): P(dart) as a function of width,
# logistic curve from the brms univariate width fit, threshold μ,
# 95% credible interval shaded around the threshold. Complete Palpa points
# are plotted along the curve at their model-predicted P(dart); incomplete
# points are excluded because their widths are lower bounds on the
# original. Hafted reference points shown along the y = 1 (darts) and
# y = 0 (arrows) edges, jittered for visibility.
# -----------------------------------------------------------------------------

grid_w_p <- data.frame(width = seq(5, 45, length.out = 1000))
grid_w_p$p_dart <- apply(posterior_epred(fit_width, newdata = grid_w_p),
                         2, median)

ref_plot <- train_lab %>%
  mutate(p_pos = ifelse(weapon == "Dart", 1.0, 0.0))

p3 <- ggplot() +
  # 95% CrI shaded around threshold
  annotate("rect",
           xmin = threshold_ci[1], xmax = threshold_ci[2],
           ymin = -Inf, ymax = Inf,
           fill = "gray70", alpha = 0.30) +
  geom_vline(xintercept = threshold_med,
             colour = "black", linetype = "dashed", linewidth = 0.5) +
  # Logistic P(dart) curve, coloured along its length with the same red-
  # yellow-blue gradient as Plot 1's background.
  geom_path(data = grid_w_p,
            aes(width, p_dart, colour = p_dart),
            linewidth = 1.2) +
  scale_colour_gradient2(
    low = pdart_low, mid = pdart_mid, high = pdart_high,
    midpoint = 0.5, limits = c(0, 1),
    name = "P(dart)", guide = "none"
  ) +
  # New colour scale for hafted-reference shape outlines
  new_scale_colour() +
  # Hafted reference at edges
  geom_point(data = ref_plot,
             aes(width, p_pos, shape = weapon),
             colour = "gray25", fill = NA,
             alpha = 0.55, size = 1.8, stroke = 0.6,
             position = position_jitter(width = 0, height = 0.04, seed = 1)) +
  scale_shape_manual(values = c("Arrow" = 1, "Dart" = 2),
                     name = "Hafted reference") +
  # Separate fill scale for Palpa points (complete only; all squares).
  geom_point(data = plot_palpa_complete,
             aes(width, P_dart_width, fill = grade_group),
             shape = 22, colour = "black", stroke = 0.5,
             size = 2.4, alpha = 0.95) +
  scale_fill_manual(
    values = palpa_fill, name = "Palpa integrity",
    guide = guide_legend(
      override.aes = list(shape = 22, colour = "black", size = 3.2)
    )
  ) +
  annotate("text",
           x = threshold_med, y = 0.52,
           label = sprintf("μ = %.2f mm\n(95%% CrI: %.2f – %.2f)",
                           threshold_med, threshold_ci[1], threshold_ci[2]),
           hjust = -0.05, vjust = 0.5, size = 3.3, colour = "black") +
  scale_x_continuous(limits = c(5, 45), breaks = seq(5, 45, 5)) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.25)) +
  labs(
    x = "Maximum width (mm)",
    y = "P(dart) — univariate width model",
    title = "Width threshold for dart/arrow classification",
    subtitle = sprintf(
      "Logistic regression on log(width); threshold μ = %.2f mm. Compare Haas and Kelly (2026): 17.6 mm.",
      threshold_med)
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("Palpa_width_logistic_plot.svg", p3,
       width = 10, height = 6, device = "svg")

cat("\n--- Plots saved ---\n")
cat("  Palpa_classification_plot.svg   (bivariate decision surface)\n")
cat("  Palpa_pdart_summary.svg         (P(dart) across the five models)\n")
cat("  Palpa_width_logistic_plot.svg   (univariate width, Haas and Kelly-comparable)\n")
cat("\nDone.\n")
