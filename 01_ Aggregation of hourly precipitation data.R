in_dir  <- "C:/Users/USER/Downloads/precipitation/data/hourly"
out_dir <- "C:/Users/USER/Downloads/precipitation/data/monthly/"

files <- list.files(in_dir, pattern = "\\.nc$", full.names = TRUE)

for (f in files) {
  
  cat("/nProcessing:", basename(f), "/n")
  
  # Load raster
  r <- rast(f)
  
  # Read time from NetCDF
  nc <- nc_open(f)
  time_raw <- nc$dim$valid_time$vals
  nc_close(nc)
  
  # Convert to POSIX time
  time <- as.POSIXct(time_raw, origin="1970-01-01", tz="UTC")
  
  # Assign names
  names(r) <- format(time, "%Y-%m-%d %H:%M")
  
  # Monthly grouping index
  month_index <- format(time, "%Y-%m")
  
  # Monthly sum
  monthly_sum <- tapp(r, month_index, fun = sum, na.rm = TRUE)
  
  names(monthly_sum) <- unique(month_index)
  
  # Build output filename
  out_file <- file.path(
    out_dir,
    paste0("monthly_", basename(f))
  )
  
  # Save output
  writeRaster(monthly_sum, out_file, overwrite = TRUE)
  
  cat("Saved:", out_file, "/n")
}