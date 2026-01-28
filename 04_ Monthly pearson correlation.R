nc_file<- "C:/Users/USER/Downloads/precipitation/data/era5_monthly_precip_m.nc"

monthly_precip<-rast(nc_file)
monthly_precip <- monthly_precip * 1000

# --- 1. Preparation and Setup ---

# Determine the total number of layers and years in your combined dataset
num_layers <- nlyr(monthly_precip)

# Note: This assumes your combined raster is perfectly aligned (e.g., 360 layers = 30 years * 12 months)
num_years <- floor(num_layers / 12) 

sea_ice <- read.csv("C:/Users/USER/Downloads/df_sea_ice.csv")

# --- 2. Correlation Loop (Month by Month) ---

# Loop iterates 12 times (for Month 1 through Month 12)
# Initialize a list to store the 12 resulting correlation maps
monthly_cor_maps <- list()
monthly_p_value_maps <- list() # New list for p-values

# Loop iterates 12 times (for Month 1 through Month 12)
for (month_num in 1:12) {
  
  # A. Identify Indices: Same as before
  monthly_indices <- seq(from = month_num, by = 12, length.out = num_years)
  
  # B. Subset Data: Same as before
  precip_subset <- monthly_precip[[monthly_indices]]
  sea_ice_subset <- sea_ice$Sea_Ice_Extent_km2[monthly_indices]
  
  # C. Prepare for Correlation (Extract Matrix)
  precip_matrix_subset <- values(precip_subset)
  
  # D. Perform Correlation and P-Value Calculation (Pixel-by-Pixel)
  
  # Initialize vectors to store the results
  correlation_vector_subset <- rep(NA, nrow(precip_matrix_subset))
  p_value_vector_subset <- rep(NA, nrow(precip_matrix_subset)) # New vector for p-values
  
  # Function to calculate both R and P in one go
  results_function <- function(x, y) {
    if (sd(x, na.rm = TRUE) == 0) {
      return(c(r = NA, p = NA))
    }
    # cor.test performs the test and returns a list object
    test_result <- cor.test(x, y, use = "complete.obs")
    # FIX: Use unname() to strip the intrinsic name "cor" from test_result$estimate.
    # This ensures the output column is correctly named 'r'
    r_value <- unname(test_result$estimate) 
    p_value <- test_result$p.value
    
    # Return a named vector with clean names 'r' and 'p'
    return(c(r = r_value, p = p_value))
  }
  
  # Apply the function across all pixels
  results_matrix <- t(apply(
    X = precip_matrix_subset, 
    MARGIN = 1, 
    FUN = function(x) results_function(x, sea_ice_subset)
  ))
  
  # Separate R and P results
  correlation_vector_subset <- results_matrix[, 'r']
  p_value_vector_subset <- results_matrix[, 'p']
  
  # E. Reconstruct the Maps and Store
  
  # 1. Correlation Map (Same as before)
  cor_map <- precip_subset[[1]]
  values(cor_map) <- correlation_vector_subset
  names(cor_map) <- paste0("Cor_", month.name[month_num])
  monthly_cor_maps[[month_num]] <- cor_map
  
  # 2. P-Value Map (NEW)
  p_map <- precip_subset[[1]]
  values(p_map) <- p_value_vector_subset
  names(p_map) <- paste0("Pval_", month.name[month_num])
  monthly_p_value_maps[[month_num]] <- p_map
}

# --- 3. Final Output (Combined Stacks) ---

# Combine the 12 correlation maps
seasonal_correlation_stack <- rast(monthly_cor_maps)

# Combine the 12 p-value maps (NEW)
seasonal_p_value_stack <- rast(monthly_p_value_maps)

# Save the final stacks
writeCDF(seasonal_correlation_stack, "seasonal_correlation_precip_seaice.nc", overwrite = TRUE)
writeCDF(seasonal_p_value_stack, "seasonal_pvalue_precip_seaice.nc", overwrite = TRUE)

print("Correlation and P-Value stacks (12 layers each) saved.")

##----------------PLOTTING--------------------------------------------------------------

# --- 1. Define Variables (Load Stacks) ---
# Assuming these files exist in the specified path:
seasonal_correlation_stack <- rast("C:/Users/USER/Downloads/precipitation/data/seasonal_correlation_precip_seaice.nc")
seasonal_p_value_stack <- rast("C:/Users/USER/Downloads/precipitation/data/seasonal_pvalue_precip_seaice.nc")

significant_points <- as.points(seasonal_p_value_stack) %>%
  st_as_sf() %>%
  st_set_geometry("geometry") %>%
  # Use the prefix found in your names() output
  pivot_longer(
    cols = starts_with("seasonal_pvalue_precip_seaice_Z1="), 
    names_to = "Month_Index",
    values_to = "p_val"
  ) %>%
  # Categorize into significance levels
  mutate(
    sig_level = case_when(
      p_val < 0.01 ~ "99%",
      p_val < 0.02 ~ "98%",
      p_val < 0.05 ~ "95%",
      TRUE ~ NA_character_
    )
  ) %>%
  # Filter out non-significant points
  filter(!is.na(sig_level)) %>%
  # Extract month number and convert to factor for correct order
  mutate(
    Month_Num = as.integer(sub(".*Z1=", "", Month_Index)),
    Month = factor(month.abb[Month_Num], levels = month.abb),
    sig_level = factor(sig_level, levels = c("95%", "98%", "99%"))
  ) %>%
  select(-Month_Index, -Month_Num)


correlation_df_long <- as.data.frame(
  seasonal_correlation_stack,
  xy = TRUE,
  na.rm = TRUE
) %>%
  pivot_longer(
    cols = starts_with("seasonal_correlation_precip_seaice_Z1="),
    names_to = "Month_Index",
    values_to = "R_Value"
  ) %>%
  mutate(
    Month_Num = as.integer(sub(".*Z1=", "", Month_Index)),
    Month = factor(month.abb[Month_Num], levels = month.abb)
  ) %>%
  select(-Month_Index, -Month_Num)

# 1b. ALL-VALUES LAYER (Original stack)
all_correlation_stack <- seasonal_correlation_stack

# --- 2. Prepare Geographical Boundaries using Local Shapefiles ---

# Load the local shapefiles (as specified in your prompt)
mys <- vect("C:/Users/USER/Downloads/gadm41_MYS_shp/gadm41_MYS_1.shp")
idn <- vect("C:/Users/USER/Downloads/gadm41_IDN_shp/gadm41_IDN_1.shp")
thd <- vect("C:/Users/USER/Downloads/gadm41_THA_shp/gadm41_THA_1.shp")
khm <- vect("C:/Users/USER/Downloads/gadm41_KHM_shp/gadm41_KHM_1.shp")

# --- 2. Correctly Combine SpatVector objects using merge() ---
# Use merge() to concatenate the two SpatVectors into a single SpatVector object.
region_sf <- dplyr::bind_rows(
  st_as_sf(mys),
  st_as_sf(idn),
  st_as_sf(thd),
  st_as_sf(khm)
)

region_sf <- st_transform(
  region_sf,
  st_crs(all_correlation_stack)
)

region_sf <- st_make_valid(region_sf)

region_sf <- st_crop(
  region_sf,
  st_as_sfc(st_bbox(all_correlation_stack))
)



# --- Corrected Step 3: Define the Cropping Polygon ---

regional_extent <- c(xmin = 84.875, xmax = 128.125, ymin = -10.125, ymax = 12.125)


# --- 3. Plotting Setup (Same as before) ---
cor_palette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

# --- 3. Plot 1: Full Correlation Map (All R-values) ---
month_groups <- list(
  "Jan_Apr" = c("Jan", "Feb", "Mar", "Apr"),
  "May_Aug" = c("May", "Jun", "Jul", "Aug"),
  "Sep_Dec" = c("Sep", "Oct", "Nov", "Dec")
)

for (group_name in names(month_groups)) {
  
  # Filter the correlation data and p-value data for current group
  current_months <- month_groups[[group_name]]
  df_sub <- correlation_df_long %>% filter(Month %in% current_months)
  
  # Prepare the p-value data for contours specifically for these months
  p_sub <- as.data.frame(seasonal_p_value_stack, xy = TRUE) %>%
    pivot_longer(cols = starts_with("seasonal_pvalue"), names_to = "Month_Idx", values_to = "p_val") %>%
    mutate(Month_Num = as.integer(sub(".*Z1=", "", Month_Idx)),
           Month = factor(month.abb[Month_Num], levels = month.abb)) %>%
    filter(Month %in% current_months)
  
  # Generate the plot
  p <- ggplot() +
    # 1. Background Heatmap
    geom_raster(data = df_sub, aes(x = x, y = y, fill = R_Value)) +
    
    # 2. Significance Contours (Using your Red, Yellow, Green logic)
    geom_contour(data = p_sub, aes(x = x, y = y, z = p_val), 
                 breaks = 0.10, color = "green", linewidth = 0.2) +
    geom_contour(data = p_sub, aes(x = x, y = y, z = p_val), 
                 breaks = 0.05, color = "yellow", linewidth = 0.3) +
    geom_contour(data = p_sub, aes(x = x, y = y, z = p_val), 
                 breaks = 0.01, color = "red", linewidth = 0.4) +
    
    # 3. Correlation Value Labels
    metR::geom_text_contour(data = df_sub, aes(x = x, y = y, z = R_Value),
                            stroke = 0.1, size = 2, color = "black", 
                            breaks = c(-0.6, -0.4, 0.4, 0.6), rotate = TRUE) +
    
    # 4. Geography and Formatting
    geom_sf(data = region_sf, fill = NA, color = "black", linewidth = 0.2, inherit.aes = FALSE) +
    scale_fill_gradientn(name = "Correlation (r)", colors = cor_palette, limits = c(-1, 1)) +
    facet_wrap(~Month, ncol = 2) + 
    coord_sf(xlim = c(84.875, 128.125), ylim = c(-10.125, 12.125), expand = FALSE) +
    
    # MAINTAINED TITLES
    labs(
      title = "seasonal-Precipitation–Sea Ice Correlation",
      subtitle = "Contours: Red=99%, Yellow=95%, Green=90% | Labels: Correlation (r)",
      x = "Longitude", y = "Latitude"
    ) +
    
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "gray80"))
  
  # Save with group-specific filename
  ggsave(paste0("Seasonal_Correlation_", group_name, ".png"), p, width = 12, height = 10, dpi = 300)
}
