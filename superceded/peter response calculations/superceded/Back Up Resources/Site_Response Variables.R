#CALCULATING THE RESPONSE VARIABLES - BIOMASS, ABUNDANCE AND FUNCTIONAL DIVERSITY
#27 sites
#USing only one dataset - EmobsCombined_Filled
#Version date: 13/8/25

# -- Load packages --
library(here)
library(tidyverse)
library(RColorBrewer)
library(hutilscpp)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(igraph)
library(reshape2)
library(rfishbase)
library(stringr)
library(forcats)
library(tibble)
library(FD)
library(stringdist)
library(cluster)  # for Gower distance if needed

# -- Load and merge your EmObs --
length_files <- list.files(here("Final Raw Data"), pattern = "3D point and length measurements.csv", 
                           full.names = TRUE)
EmobsCombined <- length_files %>% map_df(read.csv)

# -- Fix species names --
# Create column to combine genus and species names 
EmobsCombined <- EmobsCombined %>%
  mutate(across(c(Genus, Species), str_trim)) %>%
  mutate(Species_full = paste(Genus, Species, sep=" "))

# Define renaming vector
species_renaming <- c(
  "Chrysiptera brownriggi" = "Chrysiptera brownriggii",
  "Chaetodon semelon" = "Chaetodon semeion",
  "Oxycheilinus digrammus" = "Oxycheilinus digramma",
  "Heteroscarus acroptilus" = "Pomacentrus philippinus")

# Fix in EmobsCombined
EmobsCombined <- EmobsCombined %>%
  mutate(Species_full = recode(Species_full, !!!species_renaming))

any(grepl("semelon|brownriggi|digrammus|acroptilus$", EmobsCombined$Species_full))  # Should be FALSE
# Returns TRUE because 'brownriggi' is still part of 'brownriggii' – safe to ignore

# Get transect list from your dataset
transect_list <- unique(str_trim(sort(EmobsCombined$Period)))

# Get species list from your dataset
species_list <- unique(str_trim(sort(EmobsCombined$Species_full)))

# -- Filtering (RMS, X, Y and Z values) --
# There is no Z value/parameters - leaving out at the moment
# If using 5m = 5000, 7m = 7000, 8m = 8000
# Filter
EmobsFiltered <- EmobsCombined %>% filter(RMS..mm. < 20,
                                            Mid.X..mm. < 2500 & Mid.X..mm. > -2500,
                                            Mid.Y..mm. < 2500 & Mid.Y..mm. > -2500)
# Identify filtered-out rows
FilteredOut <- anti_join(EmobsCombined, EmobsFiltered)
# 13 rows are filtered out

# -- Add length to 3D points based on the median of species --
# Part 1: Calculate median lengths per species per transect
median_transect <- EmobsFiltered %>% group_by(Period, Species_full) %>%
  summarise(median_length = median(Length..mm., na.rm = TRUE), .groups = "drop")
# Code for extracting only the species with no median length calculated per transect
no_median_transect <- EmobsFiltered %>% group_by(Period, Species_full) %>%
  summarise(all_lengths_na = all(is.na(Length..mm.)), .groups = "drop") %>%
  filter(all_lengths_na)
nrow(no_median_transect)

# Code for median per site 
median_site <- EmobsFiltered %>% group_by(OpCode, Species_full) %>%
  summarise(median_length = median(Length..mm., na.rm = TRUE), .groups = "drop")

no_median_site <- EmobsFiltered %>% group_by(OpCode, Species_full) %>%
  summarise(all_lengths_na = all(is.na(Length..mm.)), .groups = "drop") %>%
  filter(all_lengths_na)
nrow(no_median_site)

# Code for median across all locations
median_location <- EmobsFiltered %>% group_by(Species_full) %>%
  summarise(median_length = median(Length..mm., na.rm = TRUE), .groups = "drop")

no_median_location <- EmobsFiltered %>% group_by(Species_full) %>%
  summarise(all_lengths_na = all(is.na(Length..mm.)), .groups = "drop") %>%
  filter(all_lengths_na)
nrow(no_median_location)

# Part 2: Add median length back to filtered working dataset
#   First - join transect-level medians
EmobsCombined_Filled <- EmobsFiltered %>%
  left_join(median_transect, by = c("Period", "Species_full")) %>%
  rename(median_transect = median_length)

#   Second - join site-level medians
EmobsCombined_Filled <- EmobsCombined_Filled %>%
  left_join(median_site, by = c("OpCode", "Species_full")) %>%
  rename(median_site = median_length)

#   Third - join location-level medians
EmobsCombined_Filled <- EmobsCombined_Filled %>%
  left_join(median_location, by = "Species_full") %>%
  rename(median_location = median_length)

#   Fourth - create a final filled length column using fallback logic
EmobsCombined_Filled <- EmobsCombined_Filled %>%
  mutate(Length_filled = case_when(
    !is.na(Length..mm.)        ~ Length..mm.,
    is.na(Length..mm.) & !is.na(median_transect) ~ median_transect,
    is.na(Length..mm.) & is.na(median_transect) & !is.na(median_site) ~ median_site,
    is.na(Length..mm.) & is.na(median_transect) & is.na(median_site) & !is.na(median_location) ~ median_location,
    TRUE ~ NA_real_))

# Record which source the filled length came from
EmobsCombined_Filled <- EmobsCombined_Filled %>%
  mutate(length_source = case_when(!is.na(Length..mm.) ~ "measured",
    is.na(Length..mm.) & !is.na(median_transect) ~ "transect",
    is.na(Length..mm.) & is.na(median_transect) & !is.na(median_site) ~ "site",
    is.na(Length..mm.) & is.na(median_transect) & is.na(median_site) & !is.na(median_location) ~ "location",
    TRUE ~ "unknown"))

# Summary count of sources
length_source_summary <- EmobsCombined_Filled %>% count(length_source) %>% arrange(desc(n))
print(length_source_summary)

# Then drop intermediate median columns - cleaner
EmobsCombined_Filled <- EmobsCombined_Filled %>% 
  select(-median_transect, -median_site, -median_location)
#Missing lengths will not be included in biomass but included in abundance
#Decided on using EmobsCombined_Filled for Abundance rather than separate

# Code to filter missing filled lengths
missing_length <- EmobsCombined_Filled %>% filter(is.na(Length_filled))
n_missing <- sum(is.na(EmobsCombined_Filled$Length_filled))
cat("Number of individuals removed due to missing length:", n_missing, "\n")
missing_length_list <- unique(missing_length$Species_full)
# There are lengths missing - n=66, 14 species - exclude
EmobsCombined_Filled <- EmobsCombined_Filled %>% filter(!is.na(Length_filled))

# Check to ensure Length_filled is all filled - prevents from proceeding if not all lengths are filled
if (any(is.na(EmobsCombined_Filled$Length_filled))) {
  stop("Some lengths remain missing after imputation. Review length_filling step.")
}

# After all lengths done, convert to cm - before calculating biomass
EmobsCombined_Filled <- EmobsCombined_Filled %>% mutate(Length_cm = Length_filled / 10)

# -- Add a and b columns --
# 1: Get a & b from FishBase
length_weight_data <- length_weight(species_list)
write.csv(length_weight_data, 
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/length_weight_fishbase.csv", 
          row.names = FALSE)

length_weight_summary <- length_weight_data %>%
  group_by(Species) %>% summarise(
    a_fishbase = mean(a, na.rm = TRUE),
    b_fishbase = mean(b, na.rm = TRUE),
    .groups = "drop")

# Join FishBase a/b to main dataset
EmobsCombined_Filled <- EmobsCombined_Filled %>%
  left_join(length_weight_summary, by = c("Species_full" = "Species"))

# 2: Add a & b from Luisa’s subset for missing species
length_weight_subset <- read.csv("length_weight_subset.csv")

# Clean duplicates in Luisa’s file (if needed - needed)
length_weight_subset_clean <- length_weight_subset %>%
  group_by(Species_full) %>% summarise(
    a_manual = mean(a, na.rm = TRUE),
    b_manual = mean(b, na.rm = TRUE),
    .groups = "drop")

# Join a_manual and b_manual
EmobsCombined_Filled <- EmobsCombined_Filled %>%
  left_join(length_weight_subset_clean, by = "Species_full")

# 3: Create initial fallback values (Diego (Manual) first, then FishBase)
EmobsCombined_Filled <- EmobsCombined_Filled %>% mutate(
  a_interim = if_else(!is.na(a_manual), a_manual, a_fishbase),
  b_interim = if_else(!is.na(b_manual), b_manual, b_fishbase))

# 4: Identify species still missing a/b
missing_ab <- EmobsCombined_Filled %>%
  filter(is.na(a_interim) | is.na(b_interim)) %>% distinct(Species_full)
print(missing_ab)
# 12 species missing a and b values

# 5: Genus-level fallback
# First extract Genus
EmobsCombined_Filled <- EmobsCombined_Filled %>% mutate(Genus = word(Species_full, 1))

# Create genus-level average table from available a_interim/b_interim
genus_avg_params <- EmobsCombined_Filled %>%
  filter(!is.na(a_interim) & !is.na(b_interim)) %>% group_by(Genus) %>% summarise(
    a_genus = mean(a_interim, na.rm = TRUE),
    b_genus = mean(b_interim, na.rm = TRUE),
    .groups = "drop")

# Join genus-level averages to dataset
EmobsCombined_Filled <- EmobsCombined_Filled %>% left_join(genus_avg_params, by = "Genus")

# Final a and b parameters use this fallback logic:
#   1. Manual parameters (a_manual, b_manual)
#   2. If missing, use FishBase (a_fishbase, b_fishbase)
#   3. If still missing, use genus-level mean (a_genus, b_genus)
#   4. If still missing, apply hard-coded overrides (e.g., from tribble)
EmobsCombined_Filled <- EmobsCombined_Filled %>% mutate(
    a_final = coalesce(a_interim, a_genus),
    b_final = coalesce(b_interim, b_genus))

# Track the source of a_final and b_final for transparency
EmobsCombined_Filled <- EmobsCombined_Filled %>%
  mutate(param_source = case_when(!is.na(a_manual) ~ "manual",
                                  is.na(a_manual) & !is.na(a_fishbase) ~ "fishbase",
                                  is.na(a_manual) & is.na(a_fishbase) & !is.na(a_genus) ~ "genus", TRUE ~ "custom"))

# Summary count of sources
param_source_summary <- EmobsCombined_Filled %>% count(param_source)
print(param_source_summary)

# 6: Final check
missing_final <- EmobsCombined_Filled %>%
  filter(is.na(a_final) | is.na(b_final)) %>% distinct(Species_full)
print(missing_final)
# 1 species detected

# Create table
manual_entries <- tribble(
  ~Species_full, ~a_manual_final, ~b_manual_final,
  "Taeniamia kagoshimana", 0.01096, 3.11)

# Join to dataset
EmobsCombined_Filled <- EmobsCombined_Filled %>%
  left_join(manual_entries, by = "Species_full") %>% mutate(
    a_final = if_else(is.na(a_final), a_manual_final, a_final),
    b_final = if_else(is.na(b_final), b_manual_final, b_final)) %>%
  select(-a_manual_final, -b_manual_final)

#  Final check
missing_after_manual <- EmobsCombined_Filled %>%
  filter(is.na(a_final) | is.na(b_final)) %>% distinct(Species_full)
print(missing_after_manual)

#  Clean up intermediate columns - cleaner
EmobsCombined_Filled <- EmobsCombined_Filled %>%
  select(-a_fishbase, -b_fishbase, -a_manual, -b_manual,
         -a_interim, -b_interim, -a_genus, -b_genus)

# Check
stopifnot(!any(is.na(EmobsCombined_Filled$a_final)))
stopifnot(!any(is.na(EmobsCombined_Filled$b_final)))

# -- Arranging sites by location --
# Get opCode list from your dataset
OpCode_list <- unique(EmobsCombined_Filled$OpCode)

# Create a lookup table - Location + Renamed sites
opcode_groups <- tribble(
  ~OpCode, ~Location, ~'Site', ~'Year of Adoption', ~'Status',
  "KAD_BIRD_1", "Kadavu", "Birdland", "2005-2009", "LMMA",    
  "KAD_CAT_1", "Kadavu", "Catherine", "2000-2004", "LMMA",  
  "KAD_COW_", "Kadavu", "NukuvouCOW", "2000-2004", "LMMA",      
  "KAD_FIGO", "Kadavu", "FijiGold", "2000-2004", "LMMA",         
  "KAD_LUPE", "Kadavu", "Lupes", "2000-2004", "LMMA", 
  "KAD_MaiDive", "Kadavu", "MaiDive", "2005-2009", "LMMA",  
  "KAD_OnoBack", "Kadavu", "Narikoso", "2005-2009", "LMMA",   
  "KAD_ROB", "Kadavu", "RobertsReef", "2000-2004", "LMMA",        
  "KAD_SO1", "Kadavu", "Soso1", "2005-2009", "LMMA",       
  "KAD_SO2", "Kadavu", "Soso2", "2005-2009", "LMMA",       
  "KAD_VURRO", "Kadavu", "Vurro", "2005-2009", "LMMA",        
  "Kua_Back", "Yasawas", "KuataBack", "Non-adopter","Qoliqoli",     
  "KUA_GOAT", "Yasawas", "KuataGoat", "Non-adopter", "Qoliqoli",        
  "KUA_MPA", "Yasawas", "KuataMPA", "Non-adopter", "Qoliqoli",        
  "KUA_noMPA", "Yasawas", "KuataNonMPA", "Non-adopter", "Qoliqoli",      
  "Lau_Totoya_BT", "Lau", "CoralGarden", "2010-2014", "Qoliqoli",    
  "Lau_Totoya_Joe", "Lau", "ReefCorner", "2010-2014", "Qoliqoli",   
  "Lau_VanuaBalavu", "Lau", "VBBay", "Non-adopter", "Qoliqoli",  
  "SAV_GNU_2009", "Savusavu", "GoldenNuggets", "2005-2009", "LMMA",     
  "SAV_SLR_1909", "Savusavu", "SplitRock", "2005-2009", "LMMA",  
  "WAI_ALI", "Savusavu", "AliceWonderland", "2010-2014", "Qoliqoli",         
  "WAI_MYS", "Savusavu", "MysteryReef", "2010-2014", "Qoliqoli", 
  "YA_COW", "Yasawas", "YaCOW", "Non-adopter", "Qoliqoli", 
  "YA_LIW", "Yasawas", "LittleWonders", "Non-adopter", "Qoliqoli", 
  "YA_MAC", "Yasawas", "MantaChannel", "Non-adopter", "Qoliqoli", 
  "YA_VAL", "Yasawas", "ValiBeach", "Non-adopter", "Qoliqoli", 
  "YA_VUI", "Yasawas", "Vuvui", "Non-adopter", "Qoliqoli")
write.csv(opcode_groups, "opcode_groups.csv", row.names = FALSE)

# -- CALCULATING BIOMASS --
# Calculate biomass for each individual
EmobsCombined_Filled <- EmobsCombined_Filled %>% mutate(biomass = a_final * Length_cm^b_final)

# Optional diagnostic checks
summary(EmobsCombined_Filled$biomass)

hist(EmobsCombined_Filled$biomass,
     main = "Histogram of Individual Fish Biomass (g)",
     xlab = "Biomass (g)", col = "skyblue", border = "white")

# Sum of biomass per transect and converted to kg/ha
EmobsSumBio <- EmobsCombined_Filled %>%
  group_by(Period, OpCode) %>%
  summarise(total_biomass = sum(biomass, na.rm = TRUE), .groups = "drop") %>%
  mutate(biomass_kg_ha = total_biomass * 0.04)
#250m² transect → 1 hectare = 10000/250 = 40
#1000 grams - 1 kg 
#Hence - 40/1000 = 0.04

EmobsSumBio <- EmobsSumBio %>% left_join(opcode_groups, by = "OpCode")
head(EmobsSumBio)
write.csv(EmobsSumBio, 
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/Sum_Biomass.csv", 
          row.names = FALSE)

# -- Arranging sites by year of adoption --
# Define the order of adoption categories (chronological)
adoption_order <- c("2000-2004", "2005-2009", "2010-2014", "Non-adopter")
# First get one row per Site with Year of Adoption (assuming it's consistent per site)
site_year <- EmobsSumBio %>%
  select(Site, `Year of Adoption`) %>%
  distinct()
# Make Year of Adoption an ordered factor
site_year$`Year of Adoption` <- factor(site_year$`Year of Adoption`, levels = adoption_order)
# Order sites by the Year of Adoption factor
site_year <- site_year %>%
  arrange(`Year of Adoption`)
# Now set Site as an ordered factor in your main data
EmobsSumBio$Site <- factor(EmobsSumBio$Site, levels = site_year$Site)
# Now plot with the reordered Site factor on x-axis
ggplot(EmobsSumBio, aes(x = Site, y = biomass_kg_ha, fill = `Year of Adoption`)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  scale_fill_viridis_d(option = "C", begin = 0, end = 0.8) +
  labs(
    title = "Fish Biomass per Site",
    x = "Site",
    y = "Biomass (kg/ha)",
    fill = "Year of Adoption") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18))

# Calculate mean and standard error for biomass by site and adoption year
biomass_stats <- EmobsSumBio %>%
  group_by(Site, `Year of Adoption`) %>%
  summarise(
    mean_biomass = mean(biomass_kg_ha, na.rm = TRUE), #Add mean biomass to DOV_Site_Properties
    se_biomass = sd(biomass_kg_ha, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot(EmobsSumBio, aes(x = Site, y = biomass_kg_ha, fill = `Year of Adoption`)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, color = "black") +
  geom_errorbar(data = biomass_stats,
                aes(y = mean_biomass, ymin = mean_biomass - se_biomass, ymax = mean_biomass + se_biomass),
                width = 0.3, color = "black", size = 0.7) +
  geom_point(data = biomass_stats, aes(y = mean_biomass), color = "red", size = 3) +
  scale_fill_viridis_d(option = "C", begin = 0, end = 0.8) +
  labs(
    title = "Fish Biomass by Site",
    subtitle = "Biomass measured in kg per hectare",
    x = "Site",
    y = "Biomass (kg/ha)",
    fill = "Year of Adoption") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13))

# -- CALCULATING ABUNDANCE --
# Calculate sum abundance per transect - units is individuals/250m2
TRANSECT_AREA_250M2 <- 250

EmobsAbundance <- EmobsCombined_Filled %>%
  group_by(Period, OpCode) %>%
  summarise(total_abundance = n(), .groups = "drop") %>%
  mutate(abundance_ind_250m2 = total_abundance)
head(EmobsAbundance)

EmobsAbundance <- EmobsAbundance %>% left_join(opcode_groups, by = "OpCode")
head(EmobsAbundance)
write.csv(EmobsAbundance, 
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/Sum_Abundance.csv", 
          row.names = FALSE)

# First get one row per Site with Year of Adoption (assuming it's consistent per site)
site_year <- EmobsAbundance %>%
  select(Site, `Year of Adoption`) %>%
  distinct()
# Make Year of Adoption an ordered factor
site_year$`Year of Adoption` <- factor(site_year$`Year of Adoption`, levels = adoption_order)
# Order sites by the Year of Adoption factor
site_year <- site_year %>%
  arrange(`Year of Adoption`)
# Now set Site as an ordered factor in your main data
EmobsAbundance$Site <- factor(EmobsAbundance$Site, levels = site_year$Site)
#Plot
ggplot(EmobsAbundance, aes(x = Site, y = abundance_ind_250m2, fill = `Year of Adoption`)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  scale_fill_viridis_d(option = "C", begin = 0, end = 0.8) +
  labs(
    title = "Fish Abundance by Site",
    x = "Site",
    y = "Abundance (individuals/250m²)",
    fill = "Year of Adoption") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18))

# Calculate mean and SE for abundance by site and year
abundance_stats <- EmobsAbundance %>%
  group_by(Site, `Year of Adoption`) %>%
  summarise(
    mean_abundance = mean(abundance_ind_250m2, na.rm = TRUE), #Add mean abundance to DOV_Site_Properties
    se_abundance = sd(abundance_ind_250m2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot(EmobsAbundance, aes(x = Site, y = abundance_ind_250m2, fill = `Year of Adoption`)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, color = "black") +
  geom_errorbar(data = abundance_stats,
                aes(y = mean_abundance, ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance),
                width = 0.3, color = "black", size = 0.7) +
  geom_point(data = abundance_stats, aes(y = mean_abundance), color = "red", size = 3) +
  scale_fill_viridis_d(option = "C", begin = 0, end = 0.8) +
  labs(
    title = "Fish Abundance by Site",
    subtitle = "Abundance measured as individuals per 250 m² transect",
    x = "Site",
    y = "Abundance (individuals/250 m²)",
    fill = "Year of Adoption") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13))

# -- Combining biomass and abundance using facetwrap --
# 1. Prepare biomass dataframe
biomass_plot <- EmobsSumBio %>%
  select(Site, `Year of Adoption`, biomass_kg_ha) %>%
  mutate(variable = "Biomass (kg/ha)",
         value = biomass_kg_ha) %>%
  select(Site, `Year of Adoption`, variable, value)

# 2. Prepare abundance dataframe
abundance_plot <- EmobsAbundance %>%
  select(Site, `Year of Adoption`, abundance_ind_250m2) %>%
  mutate(variable = "Abundance (ind/250m²)",
         value = abundance_ind_250m2) %>%
  select(Site, `Year of Adoption`, variable, value)

# 3. Combine
combined_plot_data <- bind_rows(biomass_plot, abundance_plot)

# 4. Ensure Year of Adoption is a factor with correct order
combined_plot_data$`Year of Adoption` <- factor(
  combined_plot_data$`Year of Adoption`,
  levels = c("2000-2004", "2005-2009", "2010-2014", "Non-adopter"))

# 5. Reorder Site based on Year of Adoption
site_order <- combined_plot_data %>%
  distinct(Site, `Year of Adoption`) %>%
  arrange(`Year of Adoption`) %>%
  pull(Site)

combined_plot_data$Site <- factor(combined_plot_data$Site, levels = site_order)

# 6. Plot
ggplot(combined_plot_data, aes(x = Site, y = value, fill = `Year of Adoption`)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6, color = "black") +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  scale_fill_viridis_d(option = "C", begin = 0, end = 0.8) +
  labs(
    title = "Fish Biomass and Abundance per Site by Year of Adoption",
    x = "Site",
    y = "Value",
    fill = "Year of Adoption") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom")

# -- Biomass and Abundance by Aoption --
# 1. Prepare site-level means for both metrics
site_level_means <- biomass_stats %>%
  select(Site, `Year of Adoption`, mean_biomass) %>%
  left_join(abundance_stats %>% select(Site, mean_abundance),
            by = "Site") %>%
  pivot_longer(cols = c(mean_biomass, mean_abundance),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         "mean_biomass" = "Biomass (kg/ha)",
                         "mean_abundance" = "Abundance (ind/250m²)"))

# 2. Order Year of Adoption as factor
site_level_means$`Year of Adoption` <- factor(site_level_means$`Year of Adoption`,
                                              levels = c("2000-2004", "2005-2009", "2010-2014", "Non-adopter"))

# 3. Plot with Year of Adoption on x-axis, metrics side-by-side
ggplot(site_level_means, aes(x = `Year of Adoption`, y = Value, fill = Metric)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    color = "black",  # set all dots to black
    alpha = 0.9,
    size = 2) +
  scale_fill_viridis_d(option = "C", begin = 0, end = 0.8) +
  labs(
    title = "Mean Fish Biomass and Abundance by Year of Adoption",
    x = "Year of Adoption",
    y = "Mean per Site",
    fill = "Metric") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18))

# -- Biomass and abundance by Status --
# Mean biomass per site
site_biomass <- EmobsSumBio %>%
  group_by(Site, Status) %>%
  summarise(value = mean(biomass_kg_ha, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "Biomass (kg/ha)")

# Mean abundance per site
site_abundance <- EmobsAbundance %>%
  group_by(Site, Status) %>%
  summarise(value = mean(abundance_ind_250m2, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "Abundance (ind/250m²)")

# Combine
site_summary <- bind_rows(site_biomass, site_abundance)

# Plot Status
ggplot(site_summary, aes(x = Status, y = value, fill = metric)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 2, alpha = 0.8, color = "black") +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  labs(
    title = "Mean Fish Biomass and Abundance by Status (per Site)",
    x = "Status",
    y = "Mean Value per Site",
    fill = "Metric") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(angle = 0))

# -- CALCULATING FUNCTIONAL DIVERSITY --
# Part 1: Create a species × trait matrix
# Load GASPAR trait
global_traits <- read_csv("C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/Fish_species_and_traits_gaspar.csv")
head(global_traits)

# Create species full column
global_traits <- global_traits %>% mutate(Species_full = paste(Genus, Species))

# Clean and prepare your species list
species_list_clean <- unique(str_trim(EmobsCombined_Filled$Species_full))
# From 252 to 238, 14 species dropped from EmobsFiltered to EmobsCombined_Filled

# Now match to your global trait dataset
trait_matched <- global_traits %>% filter(Species_full %in% species_list_clean)

# Check how many matched
matched_species <- unique(trait_matched$Species_full)
unmatched_species <- setdiff(species_list_clean, matched_species)
cat("Matched:", length(matched_species), "\n")
cat("Unmatched:", length(unmatched_species), "\n")
print(unmatched_species)

# Choose trait matrix
traits_selected <- trait_matched %>%
  select(Species_full, Sizec, Mobility, Activity, Schooling, Position, Diet_Mouillot_2014)

# Prepare trait matrix
traits_mat <- traits_selected %>% column_to_rownames(var = "Species_full")

# Unmatched Species
# Extract genus names from unmatched species
unmatched_genus <- str_extract(unmatched_species, "^[A-Za-z]+")

# Take mean of traits of global set
genus_trait_means <- global_traits %>%
  filter(Genus %in% unmatched_genus) %>%
  group_by(Genus) %>%
  summarise(
    Sizec = mean(Sizec, na.rm = TRUE),
    Mobility = round(mean(Mobility, na.rm = TRUE)),
    Activity = round(mean(Activity, na.rm = TRUE)),
    Schooling = round(mean(Schooling, na.rm = TRUE)),
    Position = round(mean(Position, na.rm = TRUE)),
    Diet_Mouillot_2014 = first(na.omit(Diet_Mouillot_2014)))  # take first available diet

# Rebuild unmatched with genus level
# Create a table with unmatched species and their genus
unmatched_df <- tibble(
  Species_full = unmatched_species,
  Genus = str_extract(unmatched_species, "^[A-Za-z]+"))

# Join with genus-level traits
genus_level_traits <- unmatched_df %>%
  left_join(genus_trait_means, by = "Genus") %>%
  select(-Genus)  # remove genus column if not needed

# Combine species with genus level
combined_traits <- bind_rows(
  trait_matched %>%
    select(Species_full, Sizec, Mobility, Activity, Schooling, Position, Diet_Mouillot_2014),
  genus_level_traits)

# Check
sum(is.na(combined_traits))

# Identify NAs
combined_traits %>% filter(if_any(everything(), is.na))

# Manually edit NAs
combined_traits <- combined_traits %>%
  mutate(
    Sizec = ifelse(Species_full == "Pomacentrus callainus", 9.5, Sizec),
    Schooling = ifelse(Species_full == "Pomacentrus callainus", 3, Schooling),
    Position = ifelse(Species_full == "Pomacentrus callainus", 1, Position),
    
    Sizec = ifelse(Species_full == "Pomacentrus spilotoceps", 7.8, Sizec),
    
    Sizec = ifelse(Species_full == "Cetoscarus ocellatus", 80, Sizec),
    Diet_Mouillot_2014 = ifelse(Species_full == "Cetoscarus ocellatus", "HM", Diet_Mouillot_2014),
    
#Archamia and Apogon species reference
    Sizec = ifelse(Species_full == "Taeniamia kagoshimana", 5.3, Sizec), #FB
    Mobility = ifelse(Species_full == "Taeniamia kagoshimana", 1, Mobility),
    Activity = ifelse(Species_full == "Taeniamia kagoshimana", 3, Activity),
    Schooling = ifelse(Species_full == "Taeniamia kagoshimana", 4, Schooling),
    Position = ifelse(Species_full == "Taeniamia kagoshimana", 2, Position),
    Diet_Mouillot_2014 = ifelse(Species_full == "Taeniamia kagoshimana", "PK", Diet_Mouillot_2014),)

# Check
combined_traits %>% filter(if_any(everything(), is.na))

# Edit these, just renamed
manual_trait_edits <- tribble(
  ~Species_full,                 ~Sizec, ~Mobility, ~Activity, ~Schooling, ~Position, ~Diet_Mouillot_2014,
  "Stegastes lacrymatus",        10,       1,         1,         1,         1,         "OM",
  "Plectroglyphidodon fasciolatus",      15,       1,         1,         3,        1,        "HD",
  "Zebrasoma velifer",        40,       2,         1,         2,        1,        "HD",
)

# Join manual edits to combined_traits
combined_traits <- combined_traits %>%
  full_join(manual_trait_edits, by = "Species_full", suffix = c("", ".manual")) %>%
  mutate(
    Sizec = coalesce(Sizec.manual, Sizec),
    Mobility = coalesce(Mobility.manual, Mobility),
    Activity = coalesce(Activity.manual, Activity),
    Schooling = coalesce(Schooling.manual, Schooling),
    Position = coalesce(Position.manual, Position),
    Diet_Mouillot_2014 = coalesce(Diet_Mouillot_2014.manual, Diet_Mouillot_2014)
  ) %>%
  select(Species_full, Sizec, Mobility, Activity, Schooling, Position, Diet_Mouillot_2014)

# Check
combined_traits %>% filter(if_any(-Species_full, is.na))

# Final 
traits_mat <- combined_traits %>%
  drop_na() %>%
  column_to_rownames("Species_full")

# Part 2: Create a site × species abundance matrix
# Abundance per transect x species
abundance_per_transect <- EmobsCombined_Filled %>%
  filter(!is.na(Species_full)) %>%
  group_by(OpCode, Period, Species_full) %>%
  summarise(n_individuals = n(), .groups = "drop") %>%
  mutate(abundance_250m2 = n_individuals / 1)  # 250m² per transect
head(abundance_per_transect)

# Pivot
site_species_matrix <- abundance_per_transect %>%
  pivot_wider(
    names_from = Species_full,
    values_from = abundance_250m2,
    values_fill = 0)
head(site_species_matrix)

# Convert to site × species matrix with summed biomass or abundance
site_species_matrix <- EmobsCombined_Filled %>%
  group_by(OpCode, Species_full) %>%
  summarise(Abundance = n(), .groups = "drop") %>%  # Or use sum(Biomass) if you want biomass instead
  pivot_wider(names_from = Species_full, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("OpCode")  # Makes OpCode the rownames

# Part 3: Standardize data/matrix
# 1. Get shared species
shared_species <- intersect(colnames(site_species_matrix), rownames(traits_mat))

# 2. Subset both matrices to shared species only
site_species_sub <- site_species_matrix[, shared_species]
traits_sub <- traits_mat[shared_species, ]

# Part 4: Calculate functional diversity indices using FD::dbFD()
# Convert character trait to factor
traits_sub$Diet_Mouillot_2014 <- as.factor(traits_sub$Diet_Mouillot_2014)

# Check
nrow(site_species_sub) == 27   # 27 sites
nrow(traits_sub) == ncol(site_species_sub)  # Species count matches
rownames(traits_sub) == colnames(site_species_sub)  # Species match

# Run dbFD
fd_results <- dbFD(
  x = traits_sub,
  a = site_species_sub,
  stand.x = TRUE,       # Standardize traits
  calc.FRic = TRUE,     # Ensure FRic is calculated
  calc.FDiv = TRUE,     # Ensure FDiv is calculated
  calc.CWM = TRUE,      # Include CWM
  corr = "cailliez",    # Distance correction
  m = 4,            # Use full-dimensionality - "max"
  messages = TRUE
)

# Convert to dataframe
fd_metrics <- as.data.frame(fd_results[c("FRic", "FEve", "FDiv", "FDis", "RaoQ")])
fd_metrics$Site <- rownames(fd_metrics)

# View summary
summary(fd_metrics)
write.csv(fd_metrics, 
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/Functional_Diversity_Metrics.csv", 
          row.names = FALSE)

# Merge to metadata
DOV_Sites_Properties <- read.csv("C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/DOV_Sites_Properties.csv")

# Inspect
head(fd_metrics$Site)
names(DOV_Sites_Properties)
head(DOV_Sites_Properties$Dive_Site)

# Create and match
site_lookup <- data.frame(
  Site = c("KAD_BIRD_1", "KAD_CAT_1", "KAD_COW_", "KAD_FIGO", "KAD_LUPE", 
           "KAD_MaiDive", "KAD_OnoBack", "KAD_ROB", "KAD_SO1", "KAD_SO2",       
           "KAD_VURRO", "Kua_Back", "KUA_GOAT", "KUA_MPA", "KUA_noMPA", "Lau_Totoya_BT",    
           "Lau_Totoya_Joe", "Lau_VanuaBalavu", "SAV_GNU_2009", "SAV_SLR_1909", "WAI_ALI",         
           "WAI_MYS", "YA_COW", "YA_LIW", "YA_MAC", "YA_VAL", "YA_VUI"),
  Dive_Site = c("Birdland", "Catherine", "NukuvouCOW", "FijiGold", "Lupes", 
                "Maidive", "Narikoso", "RobertsReef", "Soso1", "Soso2",
                "Vurro", "KuataBack", "KuataGoat", "KuataMPA", "KuataNonMPA", "CoralGarden", 
                "ReefCorner", "VBBay", "GoldenNuggets", "SplitRock", "AliceWonderland", 
                "MysteryReef", "YaCOW", "LittleWonders", "MantaChannel", "ValiBeach", "Vuvui"),
  Year_of_Adoption = c("2005-2009", "2000-2004", "2000-2004", "2000-2004", "2000-2004", 
                       "2005-2009", "2005-2009", "2000-2004", "2005-2009", "2005-2009", 
                       "2005-2009", "Non-adopter", "Non-adopter", "Non-adopter", "Non-adopter", "2010-2014", 
                       "2010-2014", "Non-adopter", "2005-2009", "2005-2009", "2010-2014",
                       "2010-2014", "Non-adopter", "Non-adopter", "Non-adopter", "Non-adopter", "Non-adopter")
)

# Merge
fd_metrics_full <- fd_metrics %>%
  left_join(site_lookup, by = "Site") %>%
  left_join(DOV_Sites_Properties, by = "Dive_Site")

# Check
glimpse(fd_metrics_full)
table(is.na(fd_metrics_full$lat))  # Should show 0 NAs ideally - False

# Optional - Community weighted means
fd_CWM <- fd_results$CWM
fd_CWM_df <- data.frame(OpCode = rownames(fd_CWM), fd_CWM)
fd_CWM_df <- fd_CWM_df %>%
  left_join(site_lookup, by = c("OpCode" = "Site"))
fd_CWM_df <- fd_CWM_df %>%
  left_join(DOV_Sites_Properties, by = "Dive_Site")

# Save results
write.csv(fd_metrics_full, "Site_Functional_Diversity_Metrics.csv", row.names = FALSE)
write.csv(fd_CWM_df, "Site_Community_Weighted_Means.csv", row.names = FALSE)
write.csv(EmobsCombined_Filled, "EmobsCombined_Filled.csv", row.names = FALSE)