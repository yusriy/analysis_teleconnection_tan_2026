pressure_level <- c("1000", "975", "950", "925", "900", "875", "850", 
                    "825", "800", "775", "750", "700", "650", "600", 
                    "550", "500", "450", "400", "350", "300", "250", 
                    "225", "200", "175", "150", "125", "100")

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
years <- 1990:2024

meridional_wind<- "C:/Users/USER/Downloads/precipitation/data/meridional wind"
vertical_wind <- "C:/Users/USER/Downloads/precipitation/data/vertical velocity"

file_list <- list.files(path = meridional_wind, 
                        pattern = ".*\\.nc$", 
                        full.names = TRUE)

all_meridional_wind<-list()
master_names <- c()

for (f in file_list) {
  yr <- regmatches(f, regexpr("\\d{4}", f))
  
  if (length(yr) > 0) {
    r <- rast(f)
    
    # Generate names for THIS specific year
    name_grid <- expand.grid(Level = pressure_level, Month = month_names)
    year_names <- paste0(name_grid$Month, "_", yr, "_", name_grid$Level, "hPa")
    
    # Add this year's names to the master list
    master_names <- c(master_names, year_names)
    
    all_meridional_wind[[yr]] <- r
    message("Added year: ", yr)
  }
}

if (length(all_meridional_wind) > 0) {
  
  # sort by numeric year (FIXED)
  all_meridional_wind <- all_meridional_wind[
    order(as.integer(names(all_meridional_wind)))
  ]
  
  v_combined <- rast(all_meridional_wind)
  
  # force rename (use correct object)
  stopifnot(length(master_names) == nlyr(v_combined))
  names(v_combined) <- master_names
  
  print(names(v_combined)[1:5])
}

file_list <- list.files(path = vertical_wind, 
                        pattern = ".*\\.nc$", 
                        full.names = TRUE)
all_vertical_wind<-list()
master_names <- c()

for (f in file_list) {
  yr <- regmatches(f, regexpr("\\d{4}", f))
  
  if (length(yr) > 0) {
    r <- rast(f)
    
    # Generate names for THIS specific year
    name_grid <- expand.grid(Level = pressure_level, Month = month_names)
    year_names <- paste0(name_grid$Month, "_", yr, "_", name_grid$Level, "hPa")
    
    # Add this year's names to the master list
    master_names <- c(master_names, year_names)
    
    all_vertical_wind[[yr]] <- r
    message("Added year: ", yr)
  }
}

if (length(all_vertical_wind) > 0) {
  
  # sort by numeric year (FIXED)
  all_vertical_wind <- all_vertical_wind[
    order(as.integer(names(all_vertical_wind)))
  ]
  
  vertical_wind_combined <- rast(all_vertical_wind)
  
  # force rename (use correct object)
  stopifnot(length(master_names) == nlyr(vertical_wind_combined))
  names(vertical_wind_combined) <- master_names
  
  print(names(vertical_wind_combined)[1:5])
}

v_combined
vertical_wind_combined

parse_layer_info <- function(nm_vec){
  tibble(layer = seq_along(nm_vec),
         name  = nm_vec) %>%
    mutate(
      Month = str_extract(name, "^[A-Za-z]{3}"),
      Year  = as.integer(str_extract(name, "(?<=_)[0-9]{4}(?=_)")),
      p_hPa = as.integer(str_extract(name, "(?<=_)[0-9]{2,4}(?=hPa$)"))
    )
}

v_info <- parse_layer_info(names(v_combined))
omega_info <- parse_layer_info(names(vertical_wind_combined))

sea_ice <- read.csv("C:/Users/USER/Downloads/df_sea_ice.csv")

#--------------function---------------------------------------------------------
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

build_month_section <- function(month,
                                v_combined, v_info,
                                vertical_wind_combined, omega_info,
                                high_years, low_years,
                                lon_min=90, lon_max=140, lat_min=-80, lat_max=20,
                                lat_vec = seq(-80, 20, by = 5),
                                p_vec = c(1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100),
                                vx = 0.8, k = 3.0, om_q = 0.80) {
  
  # --- v composite difference (High - Low)
  v_high <- lat_pressure_composite(v_combined, v_info, high_years, month = month,
                                   lon_min = lon_min, lon_max = lon_max,
                                   lat_min = lat_min, lat_max = lat_max)
  v_low  <- lat_pressure_composite(v_combined, v_info, low_years, month = month,
                                   lon_min = lon_min, lon_max = lon_max,
                                   lat_min = lat_min, lat_max = lat_max)
  
  v_diff <- v_high %>%
    rename(v_high = comp) %>%
    left_join(v_low %>% rename(v_low = comp), by = c("lat","p_hPa")) %>%
    mutate(v = v_high - v_low) %>%
    select(lat, p_hPa, v)
  
  # --- omega profiles for p-value + omega difference
  omega_high_y <- lat_pressure_profiles(vertical_wind_combined, omega_info, high_years,
                                        month = month, lon_min = lon_min, lon_max = lon_max,
                                        lat_min = lat_min, lat_max = lat_max)
  
  omega_low_y  <- lat_pressure_profiles(vertical_wind_combined, omega_info, low_years,
                                        month = month, lon_min = lon_min, lon_max = lon_max,
                                        lat_min = lat_min, lat_max = lat_max)
  
  omega_high_mean <- omega_high_y %>%
    group_by(lat, p_hPa) %>%
    summarise(om_high = mean(val, na.rm = TRUE), .groups = "drop")
  
  omega_low_mean <- omega_low_y %>%
    group_by(lat, p_hPa) %>%
    summarise(om_low = mean(val, na.rm = TRUE), .groups = "drop")
  
  omega_diff <- omega_high_mean %>%
    left_join(omega_low_mean, by = c("lat","p_hPa")) %>%
    mutate(omega = om_high - om_low) %>%
    select(lat, p_hPa, omega)
  
  # --- p-values (Welch t-test between groups)
  pvals_omega <- bind_rows(
    omega_high_y %>% mutate(group = "high"),
    omega_low_y  %>% mutate(group = "low")
  ) %>%
    group_by(lat, p_hPa) %>%
    summarise(
      p = tryCatch(t.test(val ~ group, var.equal = FALSE)$p.value,
                   error = function(e) NA_real_),
      .groups = "drop"
    )
  
  # --- df_plot: omega diff + p-values + v diff
  df_plot <- omega_diff %>%
    left_join(pvals_omega, by = c("lat","p_hPa")) %>%
    left_join(v_diff, by = c("lat","p_hPa"))
  
  # if missing v happens, fail early
  if (anyNA(df_plot$v)) stop("Missing v values after join for month = ", month)
  
  # --- df_grid: complete lat-p grid for filled contours
  df_grid <- df_plot %>%
    mutate(lat = round(lat, 2)) %>%
    complete(lat, p_hPa) %>%
    arrange(p_hPa, lat)
  
  # --- df_vec: vectors on a thinned grid (one arrow per lat/p cell)
  df_vec <- df_plot %>%
    mutate(lat_int = round(lat)) %>%
    filter(lat_int %in% lat_vec, p_hPa %in% p_vec) %>%
    group_by(lat_int, p_hPa) %>%
    summarise(
      lat   = mean(lat),
      p_hPa = first(p_hPa),
      v     = mean(v, na.rm = TRUE),
      omega = mean(omega, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- stable omega scaling for vector direction
  om_cap <- quantile(abs(df_vec$omega), om_q, na.rm = TRUE)
  if (!is.finite(om_cap) || om_cap == 0) om_cap <- 1e-6
  
  # choose scaling factors (tune these)
  kx <- 0.25    # degrees latitude per (m/s) after scaling
  ky <- 0.18    # vertical units in log10(p) per (scaled omega unit)
  
  df_vec <- df_vec %>%
    mutate(
      omega_cap = pmin(pmax(omega, -om_cap), om_cap),
      omega_s   = omega_cap * (vx / om_cap),   # <-- create omega_s first
      y0        = log10(p_hPa)
    ) %>%
    mutate(
      dx   = 0.50 * v,
      dy   = 0.08 * (-omega_s),
      xend = lat + dx,
      yend = y0  + dy,
      yend = pmin(pmax(yend, 2), 3)
    )
  
  list(
    month = month,
    df_plot = df_plot,
    df_grid = df_grid,
    df_vec  = df_vec
  )
}

months_vec <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
monthly_sections <- vector("list", length(months_vec))
names(monthly_sections) <- months_vec

for (m in months_vec) {
  message("Processing ", m, " ...")
  
  df_mon <- sea_ice %>% dplyr::filter(Month == m)
  low_q  <- quantile(df_mon$Sea_Ice_Extent_km2, 0.25, na.rm = TRUE)
  high_q <- quantile(df_mon$Sea_Ice_Extent_km2, 0.75, na.rm = TRUE)
  
  high_years_m <- df_mon$Year[df_mon$Sea_Ice_Extent_km2 >= high_q]
  low_years_m  <- df_mon$Year[df_mon$Sea_Ice_Extent_km2 <= low_q]
  
  monthly_sections[[m]] <- build_month_section(
    month = m,
    v_combined = v_combined, v_info = v_info,
    vertical_wind_combined = vertical_wind_combined, omega_info = omega_info,
    high_years = high_years_m,
    low_years  = low_years_m,
    lon_min = 90, lon_max = 140, lat_min = -80, lat_max = 20
  )
}


df_grid_all <- map_dfr(monthly_sections, ~ .x$df_grid,
                       .id = "month")

df_vec_all <- map_dfr(monthly_sections, ~ .x$df_vec,
                      .id = "month")
df_vec_all <- df_vec_all %>%
  mutate(
    dx = 0.50 * v,
    dy_log = -0.08 * (omega_cap / max(abs(omega_cap), na.rm = TRUE)),  # log-space scaling
    xend = lat + dx,
    pend = 10^(log10(p_hPa) + dy_log)
  )
month_levels <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

df_vec_all <- df_vec_all %>%
  mutate(month = factor(month, levels = month_levels))

df_grid_all$month <- factor(df_grid_all$month, levels = month_levels)
df_vec_all$month  <- factor(df_vec_all$month,  levels = month_levels)
lim <- max(abs(df_grid_all$omega), na.rm = TRUE)   # = 0.07480433

cb_breaks <- c(-lim, -0.05, -0.02, -0.01, 0, 0.01, 0.02, 0.05, lim)

omega_lim <- 0.08
# World map
world_map <- map_data("world")

plot_4month_panel <- function(month_subset,
                              xlim = c(-80, 20),
                              ylim = c(100, 1000),
                              bins = 15) {
  
  ggplot(
    df_grid_all %>% dplyr::filter(month %in% month_subset),
    aes(x = lat, y = p_hPa)
  ) +
    geom_contour_filled(
      aes(z = omega, fill = after_stat(level_mid)),
      bins = bins
    ) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-omega_lim, omega_lim),
      breaks = seq(-0.08, 0.08, by = 0.02),
      labels = scales::number_format(accuracy = 0.01),
      name = expression(Delta*omega~"(Pa s"^{-1}*")")
    ) +
    guides(fill = guide_colorbar(
      barheight = unit(70, "mm"),
      barwidth  = unit(5, "mm"),
      ticks = TRUE
    )) +
    
    geom_contour(aes(z = p), breaks = 0.10, colour = "green4", linewidth = 0.5) +
    geom_contour(aes(z = p), breaks = 0.05, colour = "purple", linewidth = 0.5) +
    geom_contour(aes(z = p), breaks = 0.01, colour = "red", linewidth = 0.5) +
    
    # IMPORTANT: make sure df_vec_all has y and yend (in hPa) OR change these names
    geom_segment(
      data = df_vec_all %>% dplyr::filter(month %in% month_subset),
      aes(x = lat, y = p_hPa, xend = xend, yend = pend),
      inherit.aes = FALSE,
      linewidth = 0.4,
      arrow = arrow(type = "closed", length = unit(0.04, "inches"))
    )+
    
    annotate("segment",
             x = -75, xend = -70,
             y = 120, yend = 120,
             linewidth = 0.6,
             arrow = arrow(type="closed", length=unit(0.10,"inches"))
    ) +
    annotate("text",
             x = -72.5, y = 115,
             label = "1.0 m/s",
             size = 3.2
    ) +
    
    scale_y_reverse(
      trans  = "log10",
      limits = c(100, 1000),
      breaks = c(100,200,300,500,700,850,1000),
      labels = c("100","200","300","500","700","850","1000")
    ) +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    
    scale_x_continuous(
      breaks = seq(xlim[1], xlim[2], by = 5),
      labels = function(x) {
        ifelse(x < 0, paste0(abs(x), "°S"),
               ifelse(x > 0, paste0(x, "°N"), "0°"))
      }
    ) +
    
    ggh4x::facet_wrap2(~ month, ncol = 2, axes = "all") +
    
    labs(
      x = "Latitude",
      y = "Pressure Level (hPa)",
      title = "Vertical and Meridional Circulation (80°E–120°E)",
      subtitle = expression("Shading: " * Delta * omega *
                              " (blue = enhanced ascent, red = enhanced subsidence)")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.spacing = unit(1.5, "lines"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(margin = margin(t = 3)),
      axis.title.x = element_text(margin = margin(t = 6)),
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.35),
      panel.grid.minor = element_line(colour = "grey88", linewidth = 0.25)
    )
}


p_JFMA <- plot_4month_panel(c("Jan","Feb","Mar","Apr"))
p_JFMA
p_MJJA <- plot_4month_panel(c("May","Jun","Jul","Aug"))
p_MJJA
p_SOND <- plot_4month_panel(c("Sep","Oct","Nov","Dec"))
p_SOND

ggsave("Lates Verticaland Meridional_Circulation_Jan_Apr.png", p_JFMA, width = 18, height = 16, dpi = 300)
ggsave("Latest Verticaland Meridional_Circulation_May_Aug.png", p_MJJA,  width = 18, height = 16, dpi = 300)
ggsave("Latest Verticaland Meridional_Circulation_Sep_Dec.png", p_SOND,  width = 18, height = 16, dpi = 300)

