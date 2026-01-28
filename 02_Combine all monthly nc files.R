nc_files <- list.files(
  "C:/Users/USER/Downloads/precipitation/data/monthly/",
  pattern = "\\.nc$",
  full.names = TRUE
)

nc_files

nc_files <- sort(nc_files)
tp_stack <- rast(nc_files)
writeCDF(
  tp_stack,
  "C:/Users/USER/Downloads/precipitation/data/era5_monthly_precip_m.nc",
  overwrite = TRUE
)

