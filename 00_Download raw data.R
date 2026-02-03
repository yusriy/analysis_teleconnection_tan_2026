## Open Notepad and paste this url and key 
# url: https://cds.climate.copernicus.eu/api
# key: UID:YOUR_API_KEY
# then save as file name= .cdsapirc , save as type = all files, location = "C:\Users\USER\.cdsapirc"
# API_KEY obtained from own era 5 account 
# Restart R studio
# -----------------------------------------------
# ERA5 Precipitation Download (monthly files)
# -----------------------------------------------

# 1) Install & load reticulate
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}
library(reticulate)

# 2) Make sure CDS API (Python) is available
reticulate::py_config()
reticulate::py_install("cdsapi", pip = TRUE)
cdsapi <- import("cdsapi")

# Create client
c <- cdsapi$Client()


#--------------precipitation---------------------------------------------------
for (yr in 1990:2024) {
        c$retrieve(
              "reanalysis-era5-single-levels",
              list(
                    product_type = "reanalysis",
                    variable = "total_precipitation",
                    year = as.character(yr),
                    month = sprintf("%02d", 1:12),
                    day = sprintf("%02d", 1:31),
                    time = sprintf("%02d:00", 0:23),
                    area = c(20, -43, -80, 120),
                    format = "netcdf"),
              paste0("era5_total_precip_", yr, ".nc")
    )
}

#------------------------geopotential height-----------------------------------------
# Loop through years from 1990 to 2025
for (yr in 1990:2024) {
  c$retrieve(
    "reanalysis-era5-pressure-levels-monthly-means", # Correct dataset for pressure levels
    list(
      product_type = "monthly_averaged_reanalysis",
      variable = "geopotential",
      pressure_level = list("200", "500", "850"), # Standard levels used in your provided papers
      year = as.character(yr),
      month = sprintf("%02d", 1:12),
      time = "00:00", 
      # Area [North, West, South, East] 
      # Corrected to cover India and the Indian Ocean sector of Antarctica
      area = c(20, -43, -80, 120), 
      format = "netcdf"
    ),
    paste0("era5_monthly_geopotential_", yr, ".nc")
  )
}

#--------------u,v and vertical --------------------------------------------------------
#Meridional wind
# Loop through years from 1990 to 2025
for (yr in 1990:2024) {
  c$retrieve(
    "reanalysis-era5-pressure-levels-monthly-means", 
    list(
      product_type = "monthly_averaged_reanalysis",
      variable =
        "v_component_of_wind",
      pressure_level = list("100", "125", "150",
                            "175", "200", "225",
                            "250", "300", "350",
                            "400", "450", "500",
                            "550", "600", "650",
                            "700", "750", "775",
                            "800", "825", "850",
                            "875", "900", "925",
                            "950", "975", "1000"), # Key levels for monsoon/teleconnection studies
      year = as.character(yr),
      month = sprintf("%02d", 1:12),
      time = "00:00", 
      # Corrected Area [North, West, South, East] to cover India & Southern Indian Ocean
      area = c(20, -43, -80, 120), 
      format = "netcdf"
    ),
    paste0("era5_monthly_press_meridional_wind_", yr, ".nc")
  )
}

# zonal wind
for (yr in 1990:2024) {
  c$retrieve(
    "reanalysis-era5-pressure-levels-monthly-means", 
    list(
      product_type = "monthly_averaged_reanalysis",
      variable =
        "u_component_of_wind",
      pressure_level = list("100", "125", "150",
                            "175", "200", "225",
                            "250", "300", "350",
                            "400", "450", "500",
                            "550", "600", "650",
                            "700", "750", "775",
                            "800", "825", "850",
                            "875", "900", "925",
                            "950", "975", "1000"), # Key levels for monsoon/teleconnection studies
      year = as.character(yr),
      month = sprintf("%02d", 1:12),
      time = "00:00", 
      # Corrected Area [North, West, South, East] to cover India & Southern Indian Ocean
      area = c(20, -43, -80, 120), 
      format = "netcdf"
    ),
    paste0("era5_monthly_press_zonal_wind_", yr, ".nc")
  )
}

# vertically velocity
for (yr in 1990:2024) {
  c$retrieve(
    "reanalysis-era5-pressure-levels-monthly-means", 
    list(
      product_type = "monthly_averaged_reanalysis",
      variable =
        "vertical_velocity",
      pressure_level = list("100", "125", "150",
                            "175", "200", "225",
                            "250", "300", "350",
                            "400", "450", "500",
                            "550", "600", "650",
                            "700", "750", "775",
                            "800", "825", "850",
                            "875", "900", "925",
                            "950", "975", "1000"), # Key levels for monsoon/teleconnection studies
      year = as.character(yr),
      month = sprintf("%02d", 1:12),
      time = "00:00", 
      # Corrected Area [North, West, South, East] to cover India & Southern Indian Ocean
      area = c(20, -43, -80, 120), 
      format = "netcdf"
    ),
    paste0("era5_monthly_press_vertical_velo_wind_", yr, ".nc")
  )
}

#--------------------------vertically_integrated_moisture_divergence-----------------------------------------------
# Loop through years from 1990 to 2025
for (yr in 1990:2024) {
  c$retrieve(
    "reanalysis-era5-single-levels-monthly-means", # Pre-calculated integrated data
    list(
      product_type = "monthly_averaged_reanalysis",
      variable = "vertically_integrated_moisture_divergence",
      year = as.character(yr),
      month = sprintf("%02d", 1:12),
      time = "00:00", 
      # Area [North, West, South, East] corrected for Malaysia/Indian Ocean
      area = c(40, 20, -80, 120), 
      format = "netcdf"
    ),
    paste0("era5_monthly_moisture_div_", yr, ".nc")
  )
}
