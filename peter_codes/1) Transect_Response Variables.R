#CALCULATING THE RESPONSE VARIABLES - BIOMASS, ABUNDANCE AND FUNCTIONAL DIVERSITY
#27 sites = 127 transects 
#Using only one dataset - EmobsCombined_Filled
#Transect level

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
library(dplyr)
library(ggplot2)
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
write.csv(EmobsCombined,
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/EmobsCombined.csv",
          row.names = FALSE)

# Get transect list from your dataset
transect_list <- unique(str_trim(sort(EmobsCombined$Period)))

# Get species list from your dataset
species_list <- unique(str_trim(sort(EmobsCombined$Species_full)))
write.csv(species_list,
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/EC_Species_List.csv",
          row.names = FALSE)

# Get family list from dataset
family_list <- unique(str_trim(sort(EmobsCombined$Family)))
write.csv(family_list,
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/EC_Family_List.csv",
          row.names = FALSE)

# -- Filtering (RMS, X, Y and Z values) --
# There is no Z value/parameters - leaving out at the moment
# If using 5m = 5000, 7m = 7000, 8m = 8000
# Filter
EmobsFiltered <- EmobsCombined %>% filter(RMS..mm. < 20,
                                            Mid.X..mm. < 2500 & Mid.X..mm. > -2500,
                                            Mid.Y..mm. < 2500 & Mid.Y..mm. > -2500)
write.csv(EmobsFiltered,
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/EmobsFiltered.csv",
          row.names = FALSE)

# Identify filtered-out rows
FilteredOut <- anti_join(EmobsCombined, EmobsFiltered)
# 13 rows are filtered out
write.csv(FilteredOut,
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/EmobsFiltered_List.csv",
          row.names = FALSE)

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
#Decided on using EmobsCombined_Filled for Abundance rather than separate

# Code to filter missing filled lengths
missing_length <- EmobsCombined_Filled %>% filter(is.na(Length_filled))
n_missing <- sum(is.na(EmobsCombined_Filled$Length_filled))
cat("Number of individuals removed due to missing length:", n_missing, "\n")
missing_length_list <- unique(missing_length$Species_full)
# There are lengths missing - n=66, 14 species - exclude
write.csv(missing_length_list, 
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/Excluded_Species_List.csv", 
          row.names = FALSE)
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
  ~OpCode, ~Location, ~'Site', ~'Year of Adoption', ~'Status', ~'Province', ~'Fishing Grounds',
  "KAD_BIRD_1", "Kadavu", "Birdland", "2005-2009", "LMMA", "Kadavu", "Kadavu 1",
  "KAD_CAT_1", "Kadavu", "Catherine", "2000-2004", "LMMA", "Kadavu", "Kadavu 1", 
  "KAD_COW_", "Kadavu", "NukuvouCOW", "2000-2004", "LMMA", "Kadavu", "Kadavu 2",      
  "KAD_FIGO", "Kadavu", "FijiGold", "2000-2004", "LMMA", "Kadavu", "Kadavu 1",  
  "KAD_LUPE", "Kadavu", "Lupes", "2000-2004", "LMMA", "Kadavu", "Kadavu 2", 
  "KAD_MaiDive", "Kadavu", "MaiDive", "2005-2009", "LMMA", "Kadavu", "Kadavu 1",
  "KAD_OnoBack", "Kadavu", "Narikoso", "2005-2009", "LMMA", "Kadavu", "Kadavu 1",   
  "KAD_ROB", "Kadavu", "RobertsReef", "2000-2004", "LMMA", "Kadavu", "Kadavu 2",       
  "KAD_SO1", "Kadavu", "Soso1", "2005-2009", "LMMA", "Kadavu", "Kadavu 3",   
  "KAD_SO2", "Kadavu", "Soso2", "2005-2009", "LMMA", "Kadavu", "Kadavu 3",       
  "KAD_VURRO", "Kadavu", "Vurro", "2005-2009", "LMMA", "Kadavu", "Kadavu 1",      
  "Kua_Back", "Yasawas", "KuataBack", "Non-adopter","Qoliqoli", "Ba", "Namara",     
  "KUA_GOAT", "Yasawas", "KuataGoat", "Non-adopter", "Qoliqoli", "Ba", "Namara",         
  "KUA_MPA", "Yasawas", "KuataMPA", "Non-adopter", "Qoliqoli", "Ba", "Namara",       
  "KUA_noMPA", "Yasawas", "KuataNonMPA", "Non-adopter", "Qoliqoli", "Ba", "Namara",        
  "Lau_Totoya_BT", "Lau", "CoralGarden", "2010-2014", "Qoliqoli", "Lau", "Totoya",    
  "Lau_Totoya_Joe", "Lau", "ReefCorner", "2010-2014", "Qoliqoli", "Lau", "Totoya",  
  "Lau_VanuaBalavu", "Lau", "VBBay", "Non-adopter", "Qoliqoli", "Lau", "Vanua Balavu",
  "SAV_GNU_2009", "Vanua Levu", "GoldenNuggets", "2005-2009", "LMMA", "Cakaudrove", "Savusavu",      
  "SAV_SLR_1909", "Vanua Levu", "SplitRock", "2005-2009", "LMMA", "Cakaudrove", "Savusavu",
  "WAI_ALI", "Vanua Levu", "AliceWonderland", "2010-2014", "Qoliqoli", "Cakaudrove", "Wailevu",         
  "WAI_MYS", "Vanua Levu", "MysteryReef", "2010-2014", "Qoliqoli", "Cakaudrove", "Wailevu",  
  "YA_COW", "Yasawas", "YaCOW", "Non-adopter", "Qoliqoli", "Ba", "Muaira",
  "YA_LIW", "Yasawas", "LittleWonders", "Non-adopter", "Qoliqoli", "Ba", "Muaira", 
  "YA_MAC", "Yasawas", "MantaChannel", "Non-adopter", "Qoliqoli", "Ba", "Muaira", 
  "YA_VAL", "Yasawas", "ValiBeach", "Non-adopter", "Qoliqoli", "Ba", "Muaira", 
  "YA_VUI", "Yasawas", "Vuvui", "Non-adopter", "Qoliqoli", "Ba", "Muaira")
write.csv(opcode_groups, "opcode_groups.csv", row.names = FALSE)

# Get species list from your dataset
ecf_species_list <- unique(str_trim(sort(EmobsCombined_Filled$Species_full)))
write.csv(ecf_species_list,
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/ECF_Species_List.csv",
          row.names = FALSE)

# Get family list from dataset
ecf_family_list <- unique(str_trim(sort(EmobsCombined_Filled$Family)))
write.csv(ecf_family_list,
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/ECF_Family_List.csv",
          row.names = FALSE)

# -- CALCULATING BIOMASS --
# Calculate biomass for each individual
EmobsCombined_Filled <- EmobsCombined_Filled %>% mutate(biomass_kg = a_final * Length_cm^b_final / 1000)

# Sanity check: expected units -> W(g) = a * L(cm)^b so dividing by 1000 -> kg
stopifnot(is.numeric(EmobsCombined_Filled$a_final))
# quick sanity: check a_final range and b_final range
summary(EmobsCombined_Filled$a_final)
summary(EmobsCombined_Filled$b_final)

# Optional diagnostic checks
summary(EmobsCombined_Filled$biomass_kg)

hist(EmobsCombined_Filled$biomass_kg,
     main = "Histogram of Individual Fish Biomass (g)",
     xlab = "Biomass (kg)", col = "skyblue", border = "white")

# Sum of biomass per transect and converted to kg/ha
EmobsSumBio <- EmobsCombined_Filled %>%
  group_by(Period, OpCode) %>%
  summarise(biomass_kg = sum(biomass_kg, na.rm = TRUE), .groups = "drop") %>%
  mutate(biomass_kg_ha = biomass_kg * 40) # 1000 divided by 250 = 40

EmobsSumBio <- EmobsSumBio %>% left_join(opcode_groups, by = "OpCode")
head(EmobsSumBio)
write.csv(EmobsSumBio, 
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/Sum_Biomass.csv", 
          row.names = FALSE)

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

write.csv(EmobsCombined_Filled, 
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/EmobsCombined_Filled.csv", 
          row.names = FALSE)

### ----------------------------------------------------------
### FUNCTIONAL DIVERSITY ANALYSIS — FINAL FULL SCRIPT
### ----------------------------------------------------------
### Calculates FRic, FEve, FDiv, FDis, RaoQ, and CWM
### using species × trait and transect × abundance matrices.
### Traits: GASPAR + genus-level means + manual fills.
### Combines FD results with biomass, abundance, and metadata.
### ----------------------------------------------------------
#-------------------------------------------------------------
# 1. Load and prepare trait dataset
#-------------------------------------------------------------

global_traits <- read_csv("C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/Fish_species_and_traits_gaspar.csv", show_col_types = FALSE)

global_traits <- global_traits %>%
  mutate(
    Genus = str_trim(Genus),
    Species = str_trim(Species),
    Species_full = paste(Genus, Species)
  )

#-------------------------------------------------------------
# 2. Clean and match species list to traits
#-------------------------------------------------------------

species_list_clean <- unique(str_trim(EmobsCombined_Filled$Species_full))

trait_matched <- global_traits %>%
  filter(Species_full %in% species_list_clean)

matched_species <- unique(trait_matched$Species_full)
unmatched_species <- setdiff(species_list_clean, matched_species)
cat("Matched:", length(matched_species), "\n")
cat("Unmatched:", length(unmatched_species), "\n")
print(unmatched_species)

#-------------------------------------------------------------
# 3. Genus-level trait substitution
#-------------------------------------------------------------

unmatched_genus <- str_extract(unmatched_species, "^[A-Za-z]+")

genus_trait_means <- global_traits %>%
  filter(Genus %in% unmatched_genus) %>%
  group_by(Genus) %>%
  summarise(across(c(Sizec, Mobility, Activity, Schooling, Position), ~mean(.x, na.rm = TRUE)),
            Diet_Mouillot_2014 = first(na.omit(Diet_Mouillot_2014)))

unmatched_df <- tibble(
  Species_full = unmatched_species,
  Genus = str_extract(unmatched_species, "^[A-Za-z]+")
)

genus_level_traits <- unmatched_df %>%
  left_join(genus_trait_means, by = "Genus") %>%
  mutate(trait_source = "genus") %>%
  select(-Genus)

combined_traits <- bind_rows(
  trait_matched %>%
    select(Species_full, Sizec, Mobility, Activity, Schooling, Position, Diet_Mouillot_2014),
  genus_level_traits
)

#-------------------------------------------------------------
# 4. Manual trait fixes for missing species
#-------------------------------------------------------------

combined_traits <- combined_traits %>%
  mutate(
    Sizec = ifelse(Species_full == "Pomacentrus callainus", 9.5, Sizec),
    Schooling = ifelse(Species_full == "Pomacentrus callainus", 3, Schooling),
    Position = ifelse(Species_full == "Pomacentrus callainus", 1, Position),
    
    Sizec = ifelse(Species_full == "Pomacentrus spilotoceps", 7.8, Sizec),
    
    Sizec = ifelse(Species_full == "Cetoscarus ocellatus", 80, Sizec),
    Diet_Mouillot_2014 = ifelse(Species_full == "Cetoscarus ocellatus", "HM", Diet_Mouillot_2014),
    
    # Archamia / Apogon reference species
    Sizec = ifelse(Species_full == "Taeniamia kagoshimana", 5.3, Sizec),
    Mobility = ifelse(Species_full == "Taeniamia kagoshimana", 1, Mobility),
    Activity = ifelse(Species_full == "Taeniamia kagoshimana", 3, Activity),
    Schooling = ifelse(Species_full == "Taeniamia kagoshimana", 4, Schooling),
    Position = ifelse(Species_full == "Taeniamia kagoshimana", 2, Position),
    Diet_Mouillot_2014 = ifelse(Species_full == "Taeniamia kagoshimana", "PK", Diet_Mouillot_2014)
  )

manual_trait_edits <- tribble(
  ~Species_full, ~Sizec, ~Mobility, ~Activity, ~Schooling, ~Position, ~Diet_Mouillot_2014,
  "Stegastes lacrymatus", 10, 1, 1, 1, 1, "OM",
  "Plectroglyphidodon fasciolatus", 15, 1, 1, 3, 1, "HD",
  "Zebrasoma velifer", 40, 2, 1, 2, 1, "HD"
)

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

#-------------------------------------------------------------
# 5. Prepare final trait matrix
#-------------------------------------------------------------

traits_mat <- combined_traits %>%
  drop_na() %>%
  distinct() %>%
  column_to_rownames("Species_full")

# Convert categorical traits to factors
traits_mat <- traits_mat %>%
  mutate(
    Mobility = as.factor(Mobility),
    Activity = as.factor(Activity),
    Schooling = as.factor(Schooling),
    Position = as.factor(Position),
    Diet_Mouillot_2014 = as.factor(Diet_Mouillot_2014)
  )

cat("Final trait matrix species count:", nrow(traits_mat), "\n")

#-------------------------------------------------------------
# 6. Build abundance (transect × species) matrix
#-------------------------------------------------------------

abund_matrix <- EmobsCombined_Filled %>%
  group_by(Period, Species_full) %>%
  summarise(Abundance = n(), .groups = "drop") %>%
  pivot_wider(names_from = Species_full, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("Period")

#-------------------------------------------------------------
# 7. Align trait and abundance matrices
#-------------------------------------------------------------

common_species <- sort(intersect(rownames(traits_mat), colnames(abund_matrix)))

traits_mat2 <- traits_mat[common_species, ]
abund_matrix2 <- abund_matrix[, common_species]

stopifnot(all(rownames(traits_mat2) == colnames(abund_matrix2)))

#-------------------------------------------------------------
# 8. Calculate Functional Diversity metrics (version safe)
#-------------------------------------------------------------

# Check FD version
fd_version <- as.numeric(strsplit(as.character(packageVersion("FD")), "\\.")[[1]][2])

if (fd_version < 13) {
  message("FD version ≤ 1.0-12 detected — using legacy argument set")
  fd_results <- dbFD(
    x = traits_mat2,
    a = abund_matrix2,
    stand.x = TRUE,
    corr = "cailliez",
    dist = "gower",
    calc.FRic = TRUE,
    calc.FDiv = TRUE,
    calc.CWM = TRUE,
    m = "max"
  )
} else {
  message("FD version ≥ 1.0-13 detected — using full argument set")
  fd_results <- dbFD(
    x = traits_mat2,
    a = abund_matrix2,
    stand.x = TRUE,
    corr = "cailliez",
    dist = "gower",
    calc.FRic = TRUE,
    calc.FDiv = TRUE,
    calc.FEve = TRUE,
    calc.FDis = TRUE,
    calc.RaoQ = TRUE,
    calc.CWM = TRUE,
    m = "max"
  )
}

#-------------------------------------------------------------
# 9. Extract and summarize results
#-------------------------------------------------------------

fd_metrics <- data.frame(
  Transect = rownames(abund_matrix2),
  FRic = fd_results$FRic,
  FEve = fd_results$FEve,
  FDiv = fd_results$FDiv,
  FDis = fd_results$FDis,
  RaoQ = fd_results$RaoQ
)

# Continuous trait CWMs only
cwm_df <- as.data.frame(fd_results$CWM) %>%
  select(where(is.numeric)) %>%
  rownames_to_column("Transect")

#-------------------------------------------------------------
# 10. Merge with metadata, biomass, and abundance
#-------------------------------------------------------------

DOV_Transect_Properties <- read.csv("C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/DOV_Transects_Properties.csv")

metadata_full <- DOV_Transect_Properties %>%
  left_join(EmobsSumBio %>% select(Period, biomass_kg, biomass_kg_ha),
            by = c("Transect" = "Period")) %>%
  left_join(EmobsAbundance %>% select(Period, abundance_ind_250m2),
            by = c("Transect" = "Period")) %>%
  left_join(fd_metrics, by = "Transect") %>%
  left_join(cwm_df, by = "Transect")

#-------------------------------------------------------------
# 11. Check completeness and export
#-------------------------------------------------------------

summary(metadata_full)
table(is.na(metadata_full$biomass_kg_ha))
table(is.na(metadata_full$abundance_ind_250m2))
table(is.na(metadata_full$FRic))

write.csv(metadata_full,
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/DOV_Transects_Properties_With_Biomass_Abundance_FD.csv",
          row.names = FALSE)

cat("Functional Diversity + Biomass + Abundance merged successfully.\n")