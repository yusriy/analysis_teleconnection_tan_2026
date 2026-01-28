zonal_wind<- "C:/Users/USER/Downloads/precipitation/data/zonal wind"

# Assuming you have a folder for zonal wind
file_list_u <- list.files(path = zonal_wind, pattern = ".*\\.nc$", full.names = TRUE)

# ... (Standard loading loop to create u_combined) ...

pressure_level <- c("1000", "975", "950", "925", "900", "875", "850", 
                    "825", "800", "775", "750", "700", "650", "600", 
                    "550", "500", "450", "400", "350", "300", "250", 
                    "225", "200", "175", "150", "125", "100")

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
years <- 1990:2024

all_zonal_wind<-list()
master_names <- c()

for (f in file_list_u) {
  yr <- regmatches(f, regexpr("\\d{4}", f))
  
  if (length(yr) > 0) {
    r <- rast(f)
    
    # Generate names for THIS specific year
    name_grid <- expand.grid(Level = pressure_level, Month = month_names)
    year_names <- paste0(name_grid$Month, "_", yr, "_", name_grid$Level, "hPa")
    
    # Add this year's names to the master list
    master_names <- c(master_names, year_names)
    
    all_zonal_wind[[yr]] <- r
    message("Added year: ", yr)
  }
}

if (length(all_zonal_wind) > 0) {
  
  # sort by numeric year (FIXED)
  all_zonal_wind <- all_zonal_wind[
    order(as.integer(names(all_zonal_wind)))
  ]
  
  u_combined <- rast(all_zonal_wind)
  
  # force rename (use correct object)
  stopifnot(length(master_names) == nlyr(u_combined))
  names(u_combined) <- master_names
  
  print(names(u_combined)[1:5])
}

parse_layer_info <- function(nm_vec){
  tibble(layer = seq_along(nm_vec),
         name  = nm_vec) %>%
    mutate(
      Month = str_extract(name, "^[A-Za-z]{3}"),
      Year  = as.integer(str_extract(name, "(?<=_)[0-9]{4}(?=_)")),
      p_hPa = as.integer(str_extract(name, "(?<=_)[0-9]{2,4}(?=hPa$)"))
    )
}

# Create the info table
u_info <- parse_layer_info(names(u_combined))

#---------------function-------------------------------------------
lat_pressure_composite <- function(r, info, years, month = "Feb",
                                   lon_min = 90, lon_max = 140,
                                   lat_min = -80, lat_max = 20){
  
  sel <- info %>%
    filter(Month == month, Year %in% years) %>%
    arrange(p_hPa, Year)
  
  r_sub <- r[[sel$layer]]
  r_sub <- crop(r_sub, ext(lon_min, lon_max, lat_min, lat_max))
  
  get_lat_composite <- function(layer_r){
    df <- as.data.frame(layer_r, xy = TRUE, na.rm = TRUE)
    val_col <- names(df)[3]
    df %>%
      group_by(y) %>%
      summarise(val = mean(.data[[val_col]], na.rm = TRUE), .groups = "drop") %>%
      rename(lat = y)
  }
  
  out_list <- vector("list", nlyr(r_sub))
  for(i in seq_len(nlyr(r_sub))){
    prof <- get_lat_composite(r_sub[[i]])
    out_list[[i]] <- prof %>%
      mutate(p_hPa = sel$p_hPa[i],
             Year  = sel$Year[i])
  }
  
  out <- bind_rows(out_list)
  
  # <-- THIS is what creates `comp`
  out %>%
    group_by(lat, p_hPa) %>%
    summarise(comp = mean(val, na.rm = TRUE), .groups = "drop")
}

lat_pressure_profiles <- function(r, info, years, month = "Feb",
                                  lon_min = 90, lon_max = 140,
                                  lat_min = -80, lat_max = 20){
  
  sel <- info %>%
    filter(Month == month, Year %in% years) %>%
    arrange(p_hPa, Year)
  
  r_sub <- r[[sel$layer]]
  r_sub <- crop(r_sub, ext(lon_min, lon_max, lat_min, lat_max))
  
  get_lat_profile <- function(layer_r){
    df <- as.data.frame(layer_r, xy = TRUE, na.rm = TRUE)
    val_col <- names(df)[3]
    df %>%
      group_by(y) %>%
      summarise(val = mean(.data[[val_col]], na.rm = TRUE), .groups = "drop") %>%
      rename(lat = y)
  }
  
  out_list <- vector("list", nlyr(r_sub))
  for(i in seq_len(nlyr(r_sub))){
    prof <- get_lat_profile(r_sub[[i]])
    out_list[[i]] <- prof %>%
      mutate(p_hPa = sel$p_hPa[i],
             Year  = sel$Year[i])
  }
  
  bind_rows(out_list)
}

#-------------looping----------------------------------------------------------
# List to hold monthly zonal wind results
zonal_results <- list()

for (m in month_names) {
  message("Calculating Zonal Wind Composite: ", m)
  
  # 1. Thresholds
  m_sea_ice <- sea_ice %>% filter(toupper(Month) == toupper(m))
  available_yrs <- unique(u_info$Year)
  
  h_yrs <- intersect(m_sea_ice$Year[m_sea_ice$Sea_Ice_Extent_km2 >= quantile(m_sea_ice$Sea_Ice_Extent_km2, 0.75)], available_yrs)
  l_yrs <- intersect(m_sea_ice$Year[m_sea_ice$Sea_Ice_Extent_km2 <= quantile(m_sea_ice$Sea_Ice_Extent_km2, 0.25)], available_yrs)
  
  if(length(h_yrs) < 2 | length(l_yrs) < 2) next
  
  # 2. U-Wind Composite
  u_h <- lat_pressure_composite(u_combined, u_info, h_yrs, month = m)
  u_l <- lat_pressure_composite(u_combined, u_info, l_yrs, month = m)
  u_d <- u_h %>% rename(u_h=comp) %>% 
    left_join(u_l %>% rename(u_l=comp), by=c("lat","p_hPa")) %>%
    mutate(u_anom = u_h - u_l)
  
  # 3. U-Wind Significance (Profiles for t-test)
  u_h_y <- lat_pressure_profiles(u_combined, u_info, h_yrs, month = m)
  u_l_y <- lat_pressure_profiles(u_combined, u_info, l_yrs, month = m)
  
  pvals_u <- bind_rows(u_h_y %>% mutate(g="h"), u_l_y %>% mutate(g="l")) %>%
    group_by(lat, p_hPa) %>%
    summarise(p = tryCatch(t.test(val ~ g)$p.value, error=function(e) NA), .groups="drop")
  
  # Combine
  zonal_results[[m]] <- u_d %>%
    left_join(pvals_u, by=c("lat","p_hPa")) %>%
    mutate(Month = factor(m, levels = month_names))
}

# Master dataframe for the whole year
full_u_df <- bind_rows(zonal_results)

# Define the plotting groups
u_panels <- list(
  "Jan_Apr" = c("Jan", "Feb", "Mar", "Apr"),
  "Mar_Aug" = c("May", "Jun", "Jul", "Aug"),
  "Sep_Dec" = c("Sep", "Oct", "Nov", "Dec")
)

for (p_name in names(u_panels)) {
  
  # Filter data for this panel
  panel_u <- full_u_df %>% filter(Month %in% u_panels[[p_name]])
  
  # Filter for significance dots (p < 0.05)
  stipple_u <- panel_u %>% filter(p <= 0.05, lat %% 2 == 0)
  
  p <- ggplot(panel_u, aes(x = lat, y = p_hPa)) +
    # 1. Shading: Zonal Wind Anomaly
    geom_contour_fill(aes(z = u_anom), na.fill = TRUE) +
    
    # 2. Contour Lines: Add structure to the shading
    geom_contour(aes(z = u_anom), color = "black", linewidth = 0.15) +
    
    # 3. Significance: Stippling (Black dots)
    geom_point(data = stipple_u, aes(x = lat, y = p_hPa), 
               size = 0.6, alpha = 0.4, color = "red") +
    
    # 4. Facet: 4 months per panel
    facet_wrap(~Month, ncol = 2) +
    
    # 5. Professional Meteorology Axes
    scale_y_reverse(trans = "log10", 
                    breaks = c(1000, 850, 700, 500, 300, 200, 100)) +
    scale_x_continuous(breaks = seq(-80, 60, 20), 
                       labels = function(x) paste0(abs(x), "°", ifelse(x < 0, "S", ifelse(x > 0, "N", "")))) +
    
    # 6. Color Scale (Divergent)
    scale_fill_divergent(name = expression(Delta*U~(m/s)), mid = "white") +
    
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          strip.background = element_rect(fill = "gray95")) +
    labs(title = paste("Zonal Wind (U) Composite Differences:", p_name),
         subtitle = "Shading: U-Wind Anomalies | Stippling: Significant at 95% level",
         x = "Latitude", y = "Pressure (hPa)")
  
  # Save the panel
  ggsave(paste0("Zonal_Wind_Only_", p_name, ".png"), p, width = 14, height = 10, dpi = 300)
}