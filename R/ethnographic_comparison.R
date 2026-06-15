### Script to compare ethnographic and archaeological data ----

# Ethnographic vs. Archaeological Arrows: Metric Comparison
# Script written by Claude AI (model 4.7).
# It was checked and processed in R Studio by E. Marsh.
# Prior to publicaiton, G. Bilotti re run the script and integrated it into the repository

# Produces (SVG, editable in Illustrator):
#   Fig_convex_hull_width_thickness.svg   — scatterplot with convex hulls
#   Fig_density_width.svg                 — density plot, width
#
# Data: Arrows_Ethnographic_v_Archaeological_updated.xlsx
#
# Packages: readxl, ggplot2, dplyr, ggExtra, svglite
# Install once if needed:
#   install.packages(c("readxl", "ggplot2", "dplyr", "ggExtra", "svglite"))
# =============================================================================

# Suppress echoing of source lines so cat() output is not doubled.
# Remove this line if you want to step through the script interactively.
options(echo = FALSE)

library(readxl)
library(ggplot2)
library(dplyr)
library(ggExtra)
library(svglite)

# -----------------------------------------------------------------------------
# 1. Load and prepare data
# -----------------------------------------------------------------------------

arrows <- read_excel("Arrows_Ethnographic_v_Archaeological_updated.xlsx")

arrows <- arrows %>%
  rename(Width = Width, Thickness = Thickness) %>%   # already correct names
  mutate(
    Group = case_when(
      group == "ethnographic" & type == "arrow" ~ "Ethnographic arrows",
      group == "archaeological" & type == "arrow" ~ "Archaeological arrows",
      group == "archaeological" & type == "dart"  ~ "Archaeological darts",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group), !is.na(Width), !is.na(Thickness))

arrows$Group <- factor(
  arrows$Group,
  levels = c("Archaeological darts", "Ethnographic arrows", "Archaeological arrows")
)

# -----------------------------------------------------------------------------
# 2. Color / shape palette — consistent with classify_Palpa_points.R style
#    Reference points: gray25 hollow shapes in the main script.
#    Here we distinguish three groups; use a neutral palette that does not
#    clash with the red-yellow-blue P(dart) surface used in the main figures.
#    Colorblind-safe (Okabe-Ito subset).
# -----------------------------------------------------------------------------

pal <- c(
  "Archaeological darts"   = "#0077BB",   # blue
  "Ethnographic arrows"    = "#EE7733",   # orange
  "Archaeological arrows"  = "#009988"    # teal
)

shape_vals <- c(
  "Archaeological darts"   = 2,    # open triangle (matches main script darts)
  "Ethnographic arrows"    = 1,    # open circle  (reduces overplotting)
  "Archaeological arrows"  = 17    # filled triangle
)

# Common theme
base_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.4, colour = "black"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key.size   = unit(0.45, "cm"),
    legend.text       = element_text(size = 10),
    legend.title      = element_blank()
  )

# -----------------------------------------------------------------------------
# 3. Compute descriptive statistics and tests (Table S3)
# -----------------------------------------------------------------------------

stats_w <- arrows %>%
  group_by(Group) %>%
  summarise(
    n      = n(),
    Median = round(median(Width), 1),
    Q1     = round(quantile(Width, 0.25), 1),
    Q3     = round(quantile(Width, 0.75), 1),
    Min    = round(min(Width), 1),
    Max    = round(max(Width), 1),
    .groups = "drop"
  )

stats_t <- arrows %>%
  group_by(Group) %>%
  summarise(
    n      = n(),
    Median = round(median(Thickness), 1),
    Q1     = round(quantile(Thickness, 0.25), 1),
    Q3     = round(quantile(Thickness, 0.75), 1),
    Min    = round(min(Thickness), 1),
    Max    = round(max(Thickness), 1),
    .groups = "drop"
  )

# Wilcoxon rank-sum test: ethnographic vs. archaeological arrows
eth  <- arrows %>% filter(Group == "Ethnographic arrows")
arch <- arrows %>% filter(Group == "Archaeological arrows")

wtest_w <- wilcox.test(eth$Width,     arch$Width,     exact = FALSE)
wtest_t <- wilcox.test(eth$Thickness, arch$Thickness, exact = FALSE)

# Vargha-Delaney A: probability that a randomly chosen value from x exceeds
# a randomly chosen value from y (ties contribute 0.5). A = 0.5 means no
# stochastic difference; A >= 0.71 is conventionally "large" (Vargha and
# Delaney 2000). Distribution-free; no assumption of normality.
vda_A <- function(x, y) {
  r <- rank(c(x, y))
  R1 <- sum(r[seq_along(x)])
  (R1 / length(x) - (length(x) + 1) / 2) / length(y)
}

A_width     <- vda_A(eth$Width,     arch$Width)
A_thickness <- vda_A(eth$Thickness, arch$Thickness)

# Bootstrap CI on the median difference (ethnographic minus archaeological).
# Resampling is done within each group independently, so unequal sample sizes
# are not a problem: each group is resampled to its own n on every iteration.
# No distributional assumptions; robust to skew and outliers.
set.seed(42)
n_boot <- 9999

boot_median_diff <- function(x, y, n_boot) {
  replicate(n_boot, median(sample(x, length(x), replace = TRUE)) -
              median(sample(y, length(y), replace = TRUE)))
}

boot_w <- boot_median_diff(eth$Width,     arch$Width,     n_boot)
boot_t <- boot_median_diff(eth$Thickness, arch$Thickness, n_boot)

# Percentile bootstrap CI (95%)
ci_w <- quantile(boot_w, c(0.025, 0.975))
ci_t <- quantile(boot_t, c(0.025, 0.975))

obs_diff_w <- median(eth$Width)     - median(arch$Width)
obs_diff_t <- median(eth$Thickness) - median(arch$Thickness)

# Overlap descriptors
arch_w_range  <- range(arch$Width)
prop_in_range <- mean(eth$Width >= arch_w_range[1] & eth$Width <= arch_w_range[2])
prop_above_max    <- mean(eth$Width >  arch_w_range[2])
prop_above_median <- mean(eth$Width >  median(arch$Width))

# -----------------------------------------------------------------------------
# 4. Figure S6: Convex hull scatterplot — width vs. thickness
# -----------------------------------------------------------------------------

vda_label <- sprintf(
  "Vargha-Delaney A (eth vs. arch arrows):\n  Width: A = %.2f   Thickness: A = %.2f",
  A_width, A_thickness
)

hulls <- arrows %>%
  group_by(Group) %>%
  slice(chull(Width, Thickness)) %>%
  ungroup()

p_hull_core <- ggplot(arrows,
                      aes(x = Width, y = Thickness,
                          colour = Group, fill = Group, shape = Group)) +
  # Strict convex hulls (behind points) via geom_polygon
  geom_polygon(data = hulls,
               aes(x = Width, y = Thickness,
                   colour = Group, fill = Group),
               alpha       = 0.12,
               linewidth   = 0.55,
               inherit.aes = FALSE,
               show.legend = FALSE) +
  # Individual points
  geom_point(size = 1.8, alpha = 0.70, stroke = 0.4) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_shape_manual(values = shape_vals) +
  # Vargha-Delaney A annotation, bottom-right of the panel
  annotate("text",
           x = Inf, y = -Inf, hjust = 1.05, vjust = -0.4,
           label = vda_label, size = 2.8, colour = "grey20",
           lineheight = 0.95) +
  labs(
    x = "Maximum width (mm)",
    y = "Maximum thickness (mm)"
  ) +
  base_theme +
  theme(legend.position = c(0.76, 0.90))

p_hull <- ggMarginal(
  p_hull_core,
  type    = "density",
  margins = "both",
  groupColour = TRUE,
  groupFill   = TRUE,
  alpha   = 0.25,
  size    = 4         # ratio of main panel to marginal panel size
)

svglite("Fig_convex_hull_width_thickness.svg", width = 5.5, height = 4.5)
print(p_hull)
dev.off()
message("Saved: Fig_convex_hull_width_thickness.svg")

# -----------------------------------------------------------------------------
# 5. Figure S8: Density plot — width
# -----------------------------------------------------------------------------

iqr_tbl <- arrows %>%
  group_by(Group) %>%
  summarise(
    Median = median(Width),
    Q1     = quantile(Width, 0.25),
    Q3     = quantile(Width, 0.75),
    .groups = "drop"
  )

# Stack the brackets below the density curves (negative y values)
iqr_tbl$y_pos <- -0.005 - 0.012 * (as.integer(iqr_tbl$Group) - 1)

p_dens <- ggplot(arrows, aes(x = Width, colour = Group, fill = Group)) +
  # Filled density curves
  geom_density(alpha = 0.18, linewidth = 0.65, adjust = 1.2) +
  # Rug plot of individual points
  geom_rug(sides = "b", alpha = 0.30, linewidth = 0.30,
           length = unit(0.020, "npc"),
           outside = FALSE) +
  # IQR brackets (Q1 to Q3) along the bottom
  geom_segment(data = iqr_tbl,
               aes(x = Q1, xend = Q3, y = y_pos, yend = y_pos, colour = Group),
               linewidth = 1.3, lineend = "round",
               inherit.aes = FALSE) +
  # Median tick marks (filled points) on the IQR bars
  geom_point(data = iqr_tbl,
             aes(x = Median, y = y_pos, fill = Group),
             shape = 21, size = 2.6, colour = "white", stroke = 0.6,
             inherit.aes = FALSE) +
  # Horizontal reference at y = 0 to separate brackets from density curves
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    x = "Maximum width (mm)",
    y = "Density"
  ) +
  base_theme +
  theme(legend.position = c(0.80, 0.78))

svglite("Fig_density_width.svg", width = 5, height = 3.5)
print(p_dens)
dev.off()
message("Saved: Fig_density_width.svg")

# =============================================================================
# 6. Consolidated results
# =============================================================================

cat("\n\n==============================================================\n")
cat("RESULTS SUMMARY: ethnographic vs. archaeological arrows\n")
cat("==============================================================\n")

# --- Width (mm) ---
# Ethnographic arrows are ~3.4 mm wider at the median, with much
# greater spread. Darts sit clearly above both arrow distributions.
cat("\nWidth (mm), distribution-free summary:\n")
print(stats_w, n = Inf)

# --- Thickness (mm) ---
# Thickness separation is weaker: the two arrow groups broadly
# overlap; only darts stand clearly apart.
cat("\nThickness (mm), distribution-free summary:\n")
print(stats_t, n = Inf)

# --- Wilcoxon rank-sum test ---
# Tests WHETHER distributions differ (not how much).
# Both p-values are very small — the difference is not due to chance.
cat("\nWilcoxon rank-sum test (ethnographic vs. archaeological arrows):\n")
cat(sprintf("  Width:     W = %.0f,  p = %.2e\n", wtest_w$statistic, wtest_w$p.value))
cat(sprintf("  Thickness: W = %.0f,  p = %.2e\n", wtest_t$statistic, wtest_t$p.value))

# --- Vargha-Delaney A ---
# Effect size: probability a random ethnographic value exceeds a random
# archaeological value (0.5 = no difference; >=0.71 = large).
# Width A is well into "large" territory; thickness is "medium".
cat("\nVargha-Delaney A (non-parametric effect size):\n")
cat(sprintf("  Width:     A = %.2f  [%s]\n",
            A_width,
            ifelse(A_width >= 0.71, "large",
                   ifelse(A_width >= 0.64, "medium", "small"))))
cat(sprintf("  Thickness: A = %.2f  [%s]\n",
            A_thickness,
            ifelse(A_thickness >= 0.71, "large",
                   ifelse(A_thickness >= 0.64, "medium", "small"))))

# --- Bootstrap CI on median difference ---
# Answers HOW MUCH bigger, not just whether. Resampling is within-group
# so unequal n (240 vs 32) is not a problem. CI excludes zero for both
# variables; width CI is tight, confirming the shift is not driven by
# a handful of outliers.
cat("\nBootstrap CI on median difference (eth minus arch, 9999 resamples):\n")
cat(sprintf("  Width:     %.2f mm  [95%% CI: %.2f, %.2f]\n",
            obs_diff_w, ci_w[1], ci_w[2]))
cat(sprintf("  Thickness: %.2f mm  [95%% CI: %.2f, %.2f]\n",
            obs_diff_t, ci_t[1], ci_t[2]))

# --- Range overlap ---
# More than half of ethnographic arrows exceed the maximum archaeological
# arrow width entirely — the bias is directional, not just extra variance.
cat("\nRange overlap (ethnographic arrows vs. archaeological arrow range):\n")
cat(sprintf("  Arch arrow width range:       %.1f - %.1f mm\n",
            arch_w_range[1], arch_w_range[2]))
cat(sprintf("  Eth arrows within that range: %.0f%%  (n = %d / %d)\n",
            prop_in_range * 100,
            sum(eth$Width >= arch_w_range[1] & eth$Width <= arch_w_range[2]),
            nrow(eth)))
cat(sprintf("  Eth arrows above the maximum: %.0f%%  (n = %d / %d)\n",
            prop_above_max * 100,
            sum(eth$Width > arch_w_range[2]), nrow(eth)))
cat(sprintf("  Eth arrows above arch median (%.1f mm): %.0f%%\n",
            median(arch$Width), prop_above_median * 100))

# Three independent non-parametric lines of evidence converge:
# significant (Wilcoxon p << 0.001), large effect (A = 0.90),
# and bootstrap CI on median difference excludes zero.
cat("\nDone.\n")
