# Extracting a and b values from Luisa's dataset

library(tidyverse)

length_weight_manual <- read.csv("C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/data_geb_paper_for_luisa.csv", 
                                 stringsAsFactors = FALSE)
View(length_weight_manual)

# Spilting Genus and species from underscore to space
# Create new column with formatted genus and species
length_weight_manual <- length_weight_manual %>%
  mutate(genus = sapply(strsplit(species, "_"), function(x) {
    paste0(toupper(substr(x[1], 1, 1)), substr(x[1], 2, nchar(x[1])))}),
    species_epithet = sapply(strsplit(species, "_"), function(x) x[2]))

# Combined Genus and species 
length_weight_manual <- length_weight_manual %>%
  mutate(genus_species = paste(genus, species_epithet))

#  Create a, b and species subset
length_weight_subset <- length_weight_manual %>%
  select(genus_species, coefa, coefb) %>% distinct()
head(length_weight_subset)

write.csv(length_weight_subset, 
          "C:/Users/Leslie/OneDrive/Desktop/Masters Analysis/length_weight_subset.csv", 
          row.names = FALSE)