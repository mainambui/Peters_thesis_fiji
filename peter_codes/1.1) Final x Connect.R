# Read the transect file
transect_final <- read.csv("DOV_Transects_Properties_With_Biomass_Abundance_FD.csv")
connect_data <- read.csv("DOV_connect_data.csv")
glimpse(transect_final)
glimpse(connect_data)

# Merge them: keeping all rows from DOV_site_final
transect_merged_data <- transect_final %>%
  left_join(connect_data, by = c("Dive_Site" = "site"))

# Check results
glimpse(transect_merged_data)

# (Optional) Save merged file
write.csv(transect_merged_data, "DOV_transect&connect_final.csv", row.names = FALSE)