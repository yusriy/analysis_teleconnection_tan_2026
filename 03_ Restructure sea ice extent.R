df_sea_ice <- read_excel("C:/Users/USER/Downloads/S_Sea_Ice_Index_Regional_Monthly_Data_G02135_v4.0.xlsx", sheet = "Indian-Extent-km^2",
                            skip= 1)



names(df_sea_ice)[1] <- "Year"                    # force the first column to be "Year"


#remove column start with "extent"

df_sea_ice <- df_sea_ice %>%
  
  select(-starts_with("rank"))



# 🔹 Manually rename the extent columns to months

names(df_sea_ice)[2:13] <- c("Jan_extent","Feb_extent","Mar_extent",
                                
                                "Apr_extent","May_extent","Jun_extent","Jul_extent",
                                
                                "Aug_extent","Sep_extent","Oct_extent",
                                
                                "Nov_extent","Dec_extent")



df_sea_ice <- df_sea_ice %>%
  pivot_longer( cols = -Year,              # all columns except Year
                names_to = "Month",        # create new Month column
                values_to = "Sea_Ice_Extent_km2"       # create new Extent column
                ) %>%
  mutate(Month = str_remove(Month, "_extent"),  # remove "_extent" if present
         Month = factor(Month,
                        levels = c("Jan","Feb","Mar","Apr","May","Jun",
                                   "Jul","Aug","Sep","Oct","Nov","Dec"),
                        ordered = TRUE)
         ) %>%
  arrange(Year, Month)%>%
  filter(Year>=1990,Year<=2024)

write.csv(df_sea_ice, "C:/Users/USER/Downloads/df_sea_ice.csv",
           row.names = FALSE)
