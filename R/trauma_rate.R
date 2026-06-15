### Scripts to generate trauma rates by region ----

# This script is a modification of Code S01 from Synder and Arkush (2024), available here:

# https://www.pnas.org/doi/suppl/10.1073/pnas.2410078121/suppl_file/pnas.2410078121.sd02.txt

# In order to run the code, please download and store the input CSV into the raw_data folder from here:

# https://www.pnas.org/doi/suppl/10.1073/pnas.2410078121/suppl_file/pnas.2410078121.sd01.csv

#############################################################################
# Splines per Region (N coast, C coast, S coast) + pooled spline
# Fits on logit(prop), bootstraps on logit, inverse-logits to percent
# Output: ~/Downloads/Trauma_Regions_plot.svg
# Written by ChatGPT; checked and run in R Studio by E. Marsh.
# G. Bilotti later checked the script and harmonised it for the repository.


# Load libraries
library(dplyr)
library(ggplot2)
library(scales)
library(svglite)

# helper functions
source(file.path("R", "helpers.R"))

# set up 
sourcedir <- file.path("data", "raw_data")
targetdir <- file.path("data","derived_data")

csv_file <- file.path(sourcedir, "pnas.2410078121.sd01.csv")

# define parameters (after Synder and Arkush)
regions_of_interest <- c("N coast", "C coast", "S coast")
n_boot <- 500                # bootstrap replicates (500 default)
spar_value <- 0.85
seed <- 123
min_unique_x_for_spline <- 4 # smooth.spline needs >= 4 unique x
min_rows_per_region <- 5     # require at least this many rows to attempt sensible bootstrap

set.seed(seed)

### Load data ----
if (!file.exists(csv_file)) stop("Data file not found at: ", csv_file)
trauma.dat <- read.csv(csv_file, stringsAsFactors = FALSE)

## cleaning/filtering 
required_cols <- c("Start.date","End.date","Adults.injured","N..adults.","Complexity","Region")
miss <- setdiff(required_cols, names(trauma.dat))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

numcols <- c("Start.date","End.date","Adults.injured","N..adults.")
trauma.dat[numcols] <- lapply(trauma.dat[numcols], function(x) as.numeric(as.character(x)))

# add Mid.date
trauma.dat <- trauma.dat |>
  mutate(Mid.date = (Start.date + End.date)/2)

# Smoothed proportion (k + 0.5) / (n + 1) to avoid exact 0/1
trauma.dat <- trauma.dat |>
  mutate(
    Adults.injured = ifelse(is.na(Adults.injured), NA_real_, Adults.injured),
    N..adults. = ifelse(is.na(N..adults.), NA_real_, N..adults.),
    prop_smooth = ifelse(is.na(Adults.injured) | is.na(N..adults.) | N..adults. <= 0,
                         NA_real_,
                         (Adults.injured + 0.5) / (N..adults. + 1)),
    N_adults = N..adults.
  )

# Keep only rows with Region in our interest set (but later we'll also fit pooled)
trauma.dat$Region <- as.character(trauma.dat$Region)

# Filter data for plotting (only rows with Mid.date and prop_smooth)
trauma_graph_all <- trauma.dat |>
  filter(!is.na(Mid.date), !is.na(prop_smooth), prop_smooth >= 0)

# # check if there are enough data
# if (nrow(trauma_graph_all) < 5) stop("Too few usable rows in the entire dataset (<5).")

## Compute fits for each region and pooled ----

pooled_df <- trauma_graph_all |> filter(Region %in% regions_of_interest)
# # check if there are enough data
# if (nrow(pooled_df) < 5) stop("Not enough pooled rows for a combined fit.")

# Pre-fit pooled primary to obtain pred_x grid
pooled_logit <- tryCatch({
  
  pooled_logit_vals <- qlogis(pooled_df$prop_smooth)
  
  smooth.spline(x = pooled_df$Mid.date, y = pooled_logit_vals, spar = spar_value)
  }, error = function(e) {
  
    warning("Pooled primary fit failed, will try region-based grids: ", e$message)
    NULL
    
})

pred_x_grid <- if (!is.null(pooled_logit)) pooled_logit$x else sort(unique(pooled_df$Mid.date))

# Regional results
region_results <- list()
for (r in regions_of_interest) {
  
  df_r <- trauma_graph_all |> filter(Region == r)
  
  if (nrow(df_r) < min_rows_per_region) {
    
    warning("Region '", r, "' has fewer than ", min_rows_per_region, " rows (n=", nrow(df_r), "). Skipping region fit.")
    
    region_results[[r]] <- NULL
    
    next
  }
  
  res <- fit_region_spline(df_r, pred_x = pred_x_grid, n_boot = n_boot, spar = spar_value, min_unique_x = min_unique_x_for_spline)
  
  if (is.null(res)) {
    warning("Fit failed or skipped for region: ", r)
    
    region_results[[r]] <- NULL
    
  } else {
    
    region_results[[r]] <- res |> 
      mutate(Region = r)
    
    message("Region ", r, ": rows=", nrow(df_r), "; valid_boot=", attr(res, "n_valid_boot"))
  }
}

# Combined / pooled fit (all three regions)
pooled_res <- fit_region_spline(pooled_df, pred_x = pred_x_grid, n_boot = n_boot, spar = spar_value, min_unique_x = min_unique_x_for_spline)
if (!is.null(pooled_res)) {
  pooled_res <- pooled_res |> 
    mutate(Region = "All regions")
  }

# Combine results into one data.frame for plotting
all_fits <- bind_rows(
  lapply(names(region_results), function(nm) if (!is.null(region_results[[nm]])) region_results[[nm]]),
  list(pooled_res)
)

# # check if it is empty
# if (nrow(all_fits) == 0) stop("No region fits available to plot.")

## Plot ----
# prepare data for plotting
# Filter for regions of interest
plot_points <- trauma_graph_all |>
  filter(Region %in% c(regions_of_interest, NA)) |>
  mutate(obs_pct = prop_smooth * 100)

# Colors for regions
region_colors <- c("N coast" = "midnightblue", "C coast" = "orchid1", "S coast" = "goldenrod2", "All regions" = "tomato2")
# factorise regions
all_fits$Region <- factor(all_fits$Region, levels = c(regions_of_interest, "All regions"))

p <- ggplot() +
  geom_ribbon(data = all_fits, 
              aes(x = x, ymin = ci_low, ymax = ci_high, fill = Region), alpha = 0.15) +
  geom_line(data = all_fits, 
            aes(x = x, y = fit_pct, color = Region, group = Region), linewidth = 1) +
  geom_point(data = plot_points |> filter(Region %in% regions_of_interest),
             aes(x = Mid.date, y = obs_pct, color = Region, size = pmin(N_adults, max(N_adults, na.rm = TRUE))),
             alpha = 0.8) +
  scale_color_manual(values = region_colors, name = "Region") +
  scale_fill_manual(values = region_colors, guide = "none") +
  scale_size_continuous(range = c(1, 6), name = "N (adults)") +
  scale_x_continuous(name = "Years BCE/CE", breaks = c(-3000, -2000, -1000, 0, 1000, 1550), limits = c(-3000, 1500)) +
  scale_y_continuous(name = "Percent of adults injured (%)", labels = scales::percent_format(scale = 1), limits = c(0, NA)) +
  theme_minimal(base_size = 14) +
  ggtitle("Trauma over Time by Region (spline on logit scale; percent on y-axis)") +
  theme(legend.position = "right")

# print(p)

# save plot
ggsave(file.path("figures", "Trauma_Regions_plot.svg"), p,
       width = 10, height = 8, device = "svg")