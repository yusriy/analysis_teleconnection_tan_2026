moisture_flux <-"C:/Users/USER/Downloads/precipitation/data/moisture flux"

file_list <- list.files(path = moisture_flux, 
                        pattern = ".*\\.nc$", 
                        full.names = TRUE)

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
years <- 1990:2024

all_moisture_flux<-list()
master_names <- c()

for (f in file_list) {
  yr <- regmatches(f, regexpr("\\d{4}", f))
  
  if (length(yr) > 0) {
    r <- rast(f)
    
    # Generate names for THIS specific year
    name_grid <- expand.grid( Month = month_names)
    year_names <- paste0(name_grid$Month, "_", yr)
    
    # Add this year's names to the master list
    master_names <- c(master_names, year_names)
    
    all_moisture_flux[[yr]] <- r
    message("Added year: ", yr)
  }
}

if (length(all_moisture_flux) > 0) {
  
  # sort by numeric year
  all_moisture_flux <- all_moisture_flux[
    order(as.integer(names(all_moisture_flux)))
  ]
  
  moisture_flux_combined <- rast(all_moisture_flux)
  
  # force rename (fail loudly if mismatch)
  stopifnot(length(master_names) == nlyr(moisture_flux_combined))
  names(moisture_flux_combined) <- master_names
  
  print(names(moisture_flux_combined)[1:5])
}

moisture_flux_combined

assign_time_from_depth <- function(r) {
  time_vec <- as.POSIXct(depth(r), origin = "1970-01-01", tz = "UTC")
  time(r) <- time_vec
  r
}
moisture_flux_combined <- assign_time_from_depth(moisture_flux_combined)
vimd <- moisture_flux_combined 

month_names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

#---------------function---------------------------------------------------------------
get_month_composite <- function(r, years, month_num, crop_ext = NULL) {
  r <- assign_time_from_depth(r)
  tt <- time(r)
  
  m_stack <- r[[lubridate::month(tt) == month_num]]
  m_tt <- time(m_stack)
  target <- m_stack[[lubridate::year(m_tt) %in% years]]
  
  if (!is.null(crop_ext)) target <- crop(target, crop_ext)
  
  app(target, mean, na.rm = TRUE)
}

#t-test
calc_t_test_month <- function(r, high_years, low_years, month_num, crop_ext = NULL) {
  
  r <- assign_time_from_depth(r)
  tt <- time(r)
  
  m_stack <- r[[lubridate::month(tt) == month_num]]
  m_tt <- time(m_stack)
  
  group_high <- m_stack[[lubridate::year(m_tt) %in% high_years]]
  group_low  <- m_stack[[lubridate::year(m_tt) %in% low_years]]
  
  if (!is.null(crop_ext)) {
    group_high <- crop(group_high, crop_ext)
    group_low  <- crop(group_low,  crop_ext)
  }
  
  mat_high <- values(group_high)   # [ncell x nyears_high]
  mat_low  <- values(group_low)    # [ncell x nyears_low]
  
  if (is.vector(mat_high)) mat_high <- matrix(mat_high, ncol = nlyr(group_high))
  if (is.vector(mat_low))  mat_low  <- matrix(mat_low,  ncol = nlyr(group_low))
  
  p_values <- rep(NA_real_, nrow(mat_high))
  
  for (i in seq_len(nrow(mat_high))) {
    x <- mat_high[i, ]
    y <- mat_low[i, ]
    if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2) next
    p_values[i] <- t.test(x, y)$p.value
  }
  
  p_map <- group_high[[1]]
  values(p_map) <- p_values
  p_map
}
# ========= MAIN LOOP =========
results <- vector("list", length(month_names))
names(results) <- month_names

for (m in month_names) {
  
  # 1) high/low years for this month (from your sea ice table)
  m_sea_ice <- sea_ice %>% filter(Month == m)
  
  low_thresh  <- quantile(m_sea_ice$Sea_Ice_Extent_km2, 0.25, na.rm = TRUE)
  high_thresh <- quantile(m_sea_ice$Sea_Ice_Extent_km2, 0.75, na.rm = TRUE)
  
  high_years <- m_sea_ice$Year[m_sea_ice$Sea_Ice_Extent_km2 >= high_thresh]
  low_years  <- m_sea_ice$Year[m_sea_ice$Sea_Ice_Extent_km2 <= low_thresh]
  
  # 2) month number for subsetting ERA5
  month_num <- match(m, month_names)
  
  # 3) composites + difference
  vimd_high <- get_month_composite(vimd, high_years, month_num, crop_ext = NULL)
  vimd_low  <- get_month_composite(vimd, low_years,  month_num, crop_ext = NULL)
  
  vimd_diff <- vimd_high - vimd_low
  names(vimd_diff) <- paste0("VIMD_diff_", m)
  
  df_vimd <- as.data.frame(vimd_diff, xy = TRUE, na.rm = TRUE)
  colnames(df_vimd)[3] <- "vimd_diff"
  
  # 4) p-values
  vimd_p <- calc_t_test_month(vimd, high_years, low_years, month_num, crop_ext = NULL)
  
  df_p <- as.data.frame(vimd_p, xy = TRUE, na.rm = TRUE)
  colnames(df_p)[3] <- "p_val"
  
  # 5) store outputs
  results[[m]] <- list(
    month = m,
    high_years = high_years,
    low_years = low_years,
    vimd_diff = vimd_diff,
    vimd_p = vimd_p,
    df_vimd = df_vimd,
    df_p = df_p
  )
  
  cat("Done:", m, "| High n =", length(high_years), "| Low n =", length(low_years), "\n")
}

month_names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# Combine VIMD diffs

# enforce month order for facets
month_names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

df_vimd_all <- bind_rows(lapply(month_names, function(m){
  results[[m]]$df_vimd %>% mutate(month = m)
}))

df_p_all <- bind_rows(lapply(month_names, function(m){
  results[[m]]$df_p %>% mutate(month = m)
}))

df_vimd_all$month <- factor(df_vimd_all$month, levels = month_names)
df_p_all$month    <- factor(df_p_all$month,    levels = month_names)

# keep only your map window
df_vimd_all <- df_vimd_all %>% filter(x >= 20, x <= 120, y >= -70, y <= 20)
df_p_all    <- df_p_all    %>% filter(x >= 20 ,x <= 120, y >= -70, y <= 20) 

# World map
world_map <- map_data("world")

plot_vimd_4month_panel <- function(month_subset, filename){
  
  df_v <- df_vimd_all %>% filter(month %in% month_subset)
  df_p <- df_p_all    %>% filter(month %in% month_subset)
  
  
  p <- ggplot(df_v, aes(x = x, y = y)) +
    geom_raster(aes(fill = vimd_diff)) +
    scale_fill_gradient2(
      name = "VIMD anomaly (kg m⁻²)",
      low = "blue", mid = "white", high = "red",
      midpoint = 0
    )+
    
    # IMPORTANT: do NOT remap x,y in contours; inherit from data is fine here
    geom_contour(
      data = df_p,
      aes(z = p_val),
      breaks = 0.05,
      color = "red",
      linewidth = 0.6
    )+
    
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = NA, color = "black", linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    
    
    coord_quickmap(xlim = c(20, 120), ylim = c(-70, 20), expand = FALSE) +
    
    scale_x_continuous(
      name = "Longitude",
      breaks = seq(20, 120, by = 20),
      labels = function(x) ifelse(x < 0, paste0(abs(x), "°W"),
                                  ifelse(x > 0, paste0(x, "°E"), "0°"))
    ) +
    scale_y_continuous(
      name = "Latitude",
      breaks = seq(-70, 20, by = 10),
      labels = function(y) ifelse(y < 0, paste0(abs(y), "°S"),
                                  ifelse(y > 0, paste0(y, "°N"), "0°"))
    ) +
    facet_wrap(~month, ncol = 2) +
    theme_minimal() +
    labs(
      title = " Vertical Moisture Flux Circulation: Antarctica to Malaysia (90–120°E)",
      subtitle = expression(
        "Shading: " * Delta * omega * 
          " (blue = ascent, red = subsidence) | Contours: Significance (p < 95%)"
      ))+
    theme(
      panel.spacing.x = unit(1.5, "lines"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      ),
      panel.spacing = unit(1.2, "lines"),
      strip.background = element_rect(
        fill = "grey90",
        colour = "black",
        linewidth = 0.4
      ),
      strip.text = element_text(
        face = "bold",
        size = 11
      ),
      legend.position = "right"
    )
}

periods <- list(
  Jan_Apr = c("Jan","Feb","Mar","Apr"),
  May_Aug = c("May","Jun","Jul","Aug"),
  Sep_Dec = c("Sep","Oct","Nov","Dec")
)

p1 <- plot_vimd_4month_panel(periods$Jan_Apr, "VIMD_Jan_Apr.png")
p2 <- plot_vimd_4month_panel(periods$May_Aug, "VIMD_May_Aug.png")
p3 <- plot_vimd_4month_panel(periods$Sep_Dec, "VIMD_Sep_Dec.png")

ggsave("VIMD_Jan_Apr_4panel.png", p1, width = 16, height = 12, dpi = 300)
ggsave("VIMD_May_Aug_4panel.png", p2, width = 16, height = 12, dpi = 300)
ggsave("VIMD_Sep_Dec_4panel.png", p2, width = 16, height = 12, dpi = 300)
