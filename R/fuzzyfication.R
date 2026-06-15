### Suitability script ----

# This script calculates the agricultural suitability of the study region 

# set up
if(!require("devtools")) install.packages("devtools")
devtools::install_git("https://gitlab.com/CRC1266-A2/FuzzyLandscapes.git")

library(FuzzyLandscapes)

sourcedir <- file.path("data","raw_data")
targetdir <- file.path("data","derived_data")

## load variables
# dem
dem <- terra::rast(file.path(sourcedir, "dem.tif"))
# compute slope
slope <- terra::terrain(dem, "slope")
# rivers
rivers <- sf::read_sf(file.path(sourcedir, "rivers.gpkg"), layer = "rivers_fig1") |>
  sf::st_transform(3857) |>
  terra::vect() |>
  terra::crop(dem)
# compute distance from rivers
r <- terra::rast(resolution = terra::res(dem), 
                 extent = terra::ext(dem),
                 crs = "EPSG:3857")

riv_dist <- terra::distance(
  terra::rasterize(rivers, r)
) |> terra::mask(dem)

# Save
terra::writeRaster(
  riv_dist, filename = here::here(targetdir, "rivers.tif"), 
  overwrite = TRUE)

#### Fuzzyfication ----
# Please refer to the paper and Mader et al 2026 for the full explanation of the parameters

# Terraces are generally found within 300 m from the rivers
# Sigmoid (open to one end) for river distance
fuzz_riv <- fl_create_ras(method = "sigmoid",
                          rast = riv_dist, 
                          p1 = -1,
                          p2= 300,
                          setname = "Rivers")

# elevation: most sites are between 2300 - 3500 m, with the peak being 2800 - 3200 m. 
# Use a trapezoid membership
fuzz_elev <- fl_create_ras(method = "trapezoid",
                           rast = dem,
                           p1 = 2300, 
                           p2= 2800, 
                           p3 = 3200, 
                           p4 = 3500, 
                           setname = "Elevation")

## Slope: terraces have been found (in the late intermediate period) at slopes as steep as 28 degrees
# A sigmoid is also suitable in this case
fuzz_slope <- fl_create_ras(method = "sigmoid",
                            rast = slope, 
                            p1 = -1,
                            p2= 28,
                            setname = "Slope")

## Combine the rasters
# give more weight to elevation (suitable ecozone)
suitability <- (2*fuzz_elev$FuzzyRaster + fuzz_riv$FuzzyRaster + 0.5*fuzz_slope$FuzzyRaster) -
  (2*fuzz_elev$FuzzyRaster*fuzz_riv$FuzzyRaster*(0.5*fuzz_slope$FuzzyRaster))

# export results
terra::values(suitability) <- terra::values(suitability) / max(terra::values(suitability), na.rm = TRUE)
terra::writeRaster(suitability, file.path(targetdir, "suitability.tif"), overwrite = TRUE)