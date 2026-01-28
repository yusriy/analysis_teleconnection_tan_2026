#--------------geopotential height-----------------------------------------------------
pressure_level<-c("850","500","200")
geopotential_height<-"C:/Users/USER/Downloads/precipitation/geopotential height"
all_geopotential_height<-list()
master_names <- c()
file_list <- list.files(path = geopotential_height, 
                        pattern = ".*\\.nc$", 
                        full.names = TRUE)

for (f in file_list) {
  yr <- regmatches(f, regexpr("\\d{4}", f))
  
  if (length(yr) > 0) {
    r <- rast(f)
    
    # Generate names for THIS specific year
    name_grid <- expand.grid(Level = pressure_level, Month = month_names)
    year_names <- paste0(name_grid$Month, "_", yr, "_", name_grid$Level, "hPa")
    
    # Add this year's names to the master list
    master_names <- c(master_names, year_names)
    
    all_geopotential_height[[yr]] <- r
    message("Added year: ", yr)
  }
}

# 4. Combine (Stack) all rasters and maintain names
if (length(all_geopotential_height) > 0) {
  
  # sort by numeric year
  all_geopotential_height <- all_geopotential_height[
    order(as.integer(names(all_geopotential_height)))
  ]
  
  z_combined <- rast(all_geopotential_height)
  
  # force rename (fail loudly if mismatch)
  stopifnot(length(master_names) == nlyr(z_combined))
  names(z_combined) <- master_names
  
  print(names(z_combined)[1:5])
}

z_combined

#---------------zonal&meridional wind------------------------------------------------
pressure_level <- c("1000", "975", "950", "925", "900", "875", "850", 
                    "825", "800", "775", "750", "700", "650", "600", 
                    "550", "500", "450", "400", "350", "300", "250", 
                    "225", "200", "175", "150", "125", "100")

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
years <- 1990:2024

zonal_wind<- "C:/Users/USER/Downloads/precipitation/zonal wind"
meridional_wind<-"C:/Users/USER/Downloads/precipitation/meridional wind"
all_meridional_wind <- list()
master_names <- c()
all_zonal_wind<-list()
master_names <- c()

file_list <- list.files(path = zonal_wind, 
                        pattern = ".*\\.nc$", 
                        full.names = TRUE)

for (f in file_list) {
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

# 4. Combine (Stack) all rasters and maintain names
if (length(all_zonal_wind) > 0) {
  
  # sort by numeric year
  all_zonal_wind <- all_zonal_wind[
    order(as.integer(names(all_zonal_wind)))
  ]
  
  u_combined <- rast(all_zonal_wind)
  
  # force rename (fail loudly if mismatch)
  stopifnot(length(master_names) == nlyr(u_combined))
  names(u_combined) <- master_names
  
  print(names(u_combined)[1:5])
}
v_combined
u_combined


#------------- save for 200, 500 and 850 hPa ---------------------------------------
u_200 <- subset(u_combined, grep("_200hPa", names(u_combined)))
v_200 <- subset(v_combined, grep("_200hPa", names(v_combined)))
z_200 <- subset(z_combined, grep("_200hPa", names(z_combined)))
writeCDF(z_200, "z_200_1990_2024.nc", 
         overwrite = TRUE, 
         varname = "z (Geopotential)", 
         unit = " m**2 s**-2  ")

#--------------load data--------------------------------------------------------------
## Load your sea ice data
sea_ice <- read.csv("C:/Users/USER/Downloads/df_sea_ice.csv")
# Clean month labels in CSV (prevents "all May" / missing months)
sea_ice$Month <- trimws(sea_ice$Month)
sea_ice$Month <- tools::toTitleCase(tolower(sea_ice$Month))  # "Jan"..."Dec"

# Paths (edit if needed)
z200_path <- "C:/Users/USER/Downloads/precipitation/200/z_200_1990_2024.nc"
u200_path <- "C:/Users/USER/Downloads/precipitation/200/u_200_1990_2024.nc"
v200_path <- "C:/Users/USER/Downloads/precipitation/200/v_200_1990_2024.nc"

z500_path <- "C:/Users/USER/Downloads/precipitation/500/z_500_1990_2024.nc"
u500_path <- "C:/Users/USER/Downloads/precipitation/500/u_500_1990_2024.nc"
v500_path <- "C:/Users/USER/Downloads/precipitation/500/v_500_1990_2024.nc"

z850_path <- "C:/Users/USER/Downloads/precipitation/850/z_850_1990_2024.nc"
u850_path <- "C:/Users/USER/Downloads/precipitation/850/u_850_1990_2024.nc"
v850_path <- "C:/Users/USER/Downloads/precipitation/850/v_850_1990_2024.nc"


# ------------------------------
# 1) Helper functions
# ------------------------------

# t-test p-value map for scalar raster stacks (z_high vs z_low)
calc_t_test_subset <- function(group_high, group_low) {
  mat_high <- values(group_high)
  mat_low  <- values(group_low)
  
  p_values <- sapply(seq_len(nrow(mat_high)), function(i) {
    if (sum(!is.na(mat_high[i, ])) < 2 || sum(!is.na(mat_low[i, ])) < 2) return(NA_real_)
    t.test(mat_high[i, ], mat_low[i, ])$p.value
  })
  
  p_map <- group_high[[1]]
  values(p_map) <- p_values
  p_map
}

# t-test p-value map for wind significance, using wind SPEED
calc_t_test_speed <- function(u_high, v_high, u_low, v_low) {
  sp_high <- sqrt(u_high^2 + v_high^2)
  sp_low  <- sqrt(u_low^2  + v_low^2)
  
  mat_high <- values(sp_high)
  mat_low  <- values(sp_low)
  
  p_values <- sapply(seq_len(nrow(mat_high)), function(i) {
    if (sum(!is.na(mat_high[i, ])) < 2 || sum(!is.na(mat_low[i, ])) < 2) return(NA_real_)
    t.test(mat_high[i, ], mat_low[i, ])$p.value
  })
  
  p_map <- sp_high[[1]]
  values(p_map) <- p_values
  p_map
}

# Convert NetCDF time to POSIXct and assign to raster
assign_time_from_depth <- function(z_rast, u_rast, v_rast, origin = "1970-01-01", tz = "UTC") {
  time_vec <- as.POSIXct(depth(v_rast), origin = origin, tz = tz)
  time(z_rast) <- time_vec
  time(u_rast) <- time_vec
  time(v_rast) <- time_vec
  list(z = z_rast, u = u_rast, v = v_rast)
}

# Process one pressure level -> returns monthly df (12 months)
process_level_monthly <- function(z_path, u_path, v_path,
                                  sea_ice_df,
                                  level_label,
                                  origin_time = "1970-01-01") {
  
  message("Reading level: ", level_label)
  
  z_rast <- rast(z_path)
  u_rast <- rast(u_path)
  v_rast <- rast(v_path)
  
  # Convert geopotential to geopotential height (m): z (m^2/s^2) / g0
  g0 <- 9.80665
  z_rast <- z_rast / g0
  
  # Assign time
  tmp <- assign_time_from_depth(z_rast, u_rast, v_rast, origin = origin_time, tz = "UTC")
  z_rast <- tmp$z; u_rast <- tmp$u; v_rast <- tmp$v
  
  all_months_results <- lapply(1:12, function(m) {
    
    mon_name <- month.abb[m]
    message("  Processing month: ", mon_name)
    
    # Sea ice selection for this month
    df_mon <- sea_ice_df %>% dplyr::filter(Month == mon_name)
    if (nrow(df_mon) == 0) return(NULL)
    
    low_q  <- quantile(df_mon$Sea_Ice_Extent_km2, 0.25, na.rm = TRUE)
    high_q <- quantile(df_mon$Sea_Ice_Extent_km2, 0.75, na.rm = TRUE)
    h_yrs  <- df_mon$Year[df_mon$Sea_Ice_Extent_km2 >= high_q]
    l_yrs  <- df_mon$Year[df_mon$Sea_Ice_Extent_km2 <= low_q]
    
    # Month subset in ERA5
    idx_mon <- which(lubridate::month(time(z_rast)) == m)
    if (length(idx_mon) == 0) return(NULL)
    
    z_mon_all <- z_rast[[idx_mon]]
    u_mon_all <- u_rast[[idx_mon]]
    v_mon_all <- v_rast[[idx_mon]]
    
    # Year subset for high/low groups
    yrs_mon <- lubridate::year(time(z_mon_all))
    idx_high <- which(yrs_mon %in% h_yrs)
    idx_low  <- which(yrs_mon %in% l_yrs)
    
    if (length(idx_high) < 2 || length(idx_low) < 2) return(NULL)
    
    z_high <- z_mon_all[[idx_high]]
    z_low  <- z_mon_all[[idx_low]]
    
    # Composite differences
    z_diff <- mean(z_high, na.rm = TRUE) - mean(z_low, na.rm = TRUE)
    u_diff <- mean(u_mon_all[[idx_high]], na.rm = TRUE) - mean(u_mon_all[[idx_low]], na.rm = TRUE)
    v_diff <- mean(v_mon_all[[idx_high]], na.rm = TRUE) - mean(v_mon_all[[idx_low]], na.rm = TRUE)
    
    # Significance (height + wind speed)
    z_p_val <- calc_t_test_subset(z_high, z_low)
    
    u_high <- u_mon_all[[idx_high]]; u_low <- u_mon_all[[idx_low]]
    v_high <- v_mon_all[[idx_high]]; v_low <- v_mon_all[[idx_low]]
    w_p_val <- calc_t_test_speed(u_high, v_high, u_low, v_low)
    
    # To data frames
    df_z  <- as.data.frame(z_diff, xy = TRUE) %>% dplyr::rename(z_anom = 3)
    df_w  <- as.data.frame(c(u_diff, v_diff), xy = TRUE) %>% dplyr::rename(u = 3, v = 4)
    df_p  <- as.data.frame(z_p_val, xy = TRUE) %>% dplyr::rename(p_val = 3)
    df_wp <- as.data.frame(w_p_val, xy = TRUE) %>% dplyr::rename(p_wind = 3)
    
    # IMPORTANT: prevent join mismatches due to floating coords
    df_z  <- df_z  %>% mutate(x = round(x, 4), y = round(y, 4))
    df_w  <- df_w  %>% mutate(x = round(x, 4), y = round(y, 4))
    df_p  <- df_p  %>% mutate(x = round(x, 4), y = round(y, 4))
    df_wp <- df_wp %>% mutate(x = round(x, 4), y = round(y, 4))
    
    df_z %>%
      left_join(df_w,  by = c("x", "y")) %>%
      left_join(df_p,  by = c("x", "y")) %>%
      left_join(df_wp, by = c("x", "y")) %>%
      mutate(
        month = factor(mon_name, levels = month.abb),
        level = factor(level_label, levels = c("200", "500", "850"))
      )
  })
  
  bind_rows(all_months_results)
}

# ------------------------------
# 2) Build datasets (ALL 12 months for each level)
# ------------------------------
df200 <- process_level_monthly(z200_path, u200_path, v200_path, sea_ice, level_label = "200")
df500 <- process_level_monthly(z500_path, u500_path, v500_path, sea_ice, level_label = "500")
df850 <- process_level_monthly(z850_path, u850_path, v850_path, sea_ice, level_label = "850")


df_all_levels <- bind_rows(df200, df500, df850)

# ------------------------------
# 3) Plotting: 4 months in one panel (per level)
# ------------------------------
# World map
world_map <- map_data("world")

periods <- list(
  Jan_Apr = month.abb[1:4],
  May_Aug = month.abb[5:8],
  Sep_Dec = month.abb[9:12]
)

make_4month_panel <- function(data, level_pick, month_list, title_text,
                              zlim = 25,
                              sig_level = 0.05,
                              wind_step_deg = 4,
                              wind_scale = 0.05,
                              xlim = c(-43, 120),
                              ylim = c(-70, 20)) {
  
  df_period <- data %>%
    filter(level == level_pick, month %in% month_list)
  
  if (nrow(df_period) == 0) {
    stop("No data found for level=", level_pick, " and months=", paste(month_list, collapse = ", "))
  }
  
  # Clean significance contour from binary mask (p < sig_level)
  df_sigmask <- df_period %>%
    mutate(sig = if_else(!is.na(p_val) & p_val < sig_level, 1, 0),
           x_bin = round(x / wind_step_deg) * wind_step_deg,
           y_bin = round(y / wind_step_deg) * wind_step_deg) %>%
    group_by(month, x_bin, y_bin) %>%
    summarise(x = mean(x), y = mean(y),
              sig = mean(sig, na.rm = TRUE),
              .groups = "drop")
  
  # Thin wind vectors + mark significant winds (p_wind)
  df_wind <- df_period %>%
    mutate(
      x_bin = round(x / wind_step_deg) * wind_step_deg,
      y_bin = round(y / wind_step_deg) * wind_step_deg
    ) %>%
    group_by(month, x_bin, y_bin) %>%
    summarise(
      x = mean(x), y = mean(y),
      u = mean(u, na.rm = TRUE),
      v = mean(v, na.rm = TRUE),
      p_wind = mean(p_wind, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(sig_wind = !is.na(p_wind) & p_wind < sig_level)
  
  ggplot(df_period, aes(x = x, y = y)) +
    geom_raster(aes(fill = z_anom)) +
    scale_fill_gradient2(
      name = "Geopotential anomaly (m)",
      low = "blue", mid = "white", high = "red",
      midpoint = 0
    ) +
    
    # Height significance contour
    geom_contour(
      data = df_sigmask,
      aes(x = x, y = y, z = sig),
      breaks = 0.5,
      color = "red",
      linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    
    # Wind vectors (scaled for map) — gray all + black significant
    geom_vector(
      data = df_wind,
      aes(x = x, y = y, dx = u * wind_scale, dy = v * wind_scale),
      inherit.aes = FALSE,
      color = "black",
      size = 0.4,
      arrow.length = 0.2
    ) +
    
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = NA, color = "black", linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    scale_x_continuous(
      name = "Longitude",
      limits = c(-50, 120),
      breaks = seq(-50, 120, by = 20),
      labels = function(x) {
        ifelse(x < 0, paste0(abs(x), "°W"),
               ifelse(x > 0, paste0(x, "°E"), "0°"))
      }
    ) +
    scale_y_continuous(
      name = "Latitude",
      limits = c(-70, 20),
      breaks = seq(-70, 20, by = 10),
      labels = function(y) {
        ifelse(y < 0, paste0(abs(y), "°S"),
               ifelse(y > 0, paste0(y, "°N"), "0°"))
      }
    ) +
    facet_wrap(~month, ncol = 2) +
    labs(
      title = title_text,
      subtitle = paste0("Level: ", level_pick, " hPa | Contour: p < ", sig_level),
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      panel.spacing.x = unit(1.5, "lines"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "right"
    )
}

#------------- change the pressure level 200,500,850hPa---------------------------------------------------------------------------------------
p200_JA <- make_4month_panel(df_all_levels, "200", periods$Jan_Apr,
                             "Geopotential Height Composite Difference (High SIE – Low SIE): Jan–Apr",
                             wind_step_deg = 4, wind_scale = 0.06)

p200_MA <- make_4month_panel(df_all_levels, "200", periods$May_Aug,
                             "Geopotential Height Composite Difference (High SIE – Low SIE): May-Aug",
                             wind_step_deg = 4, wind_scale = 0.06)

p200_JA <- make_4month_panel(df_all_levels, "200", periods$Jan_Apr,
                             "Geopotential Height Composite Difference (High SIE – Low SIE): Jan–Apr",
                             wind_step_deg = 4, wind_scale = 0.06)

ggsave("200hPa Geo & Hori_Jan_Apr.png", p200_JA, width = 16, height = 12, dpi = 300)
ggsave("200hPa Geo & Hori_May_Aug.png", p200_MA, width = 16, height = 12, dpi = 300)
ggsave("200hPa Geo & Hori_Sep_Dec.png", p200_SD, width = 16, height = 12, dpi = 300)

