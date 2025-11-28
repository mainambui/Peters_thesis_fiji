library(tidyverse)
library(RColorBrewer)
library(viridis)
library(patchwork)

# Read file 
dov_final <- read.csv("DOV_transect&connect_final.csv")

# Prepare data in long format
dov_long <- dov_final %>%
  select(Fishing.Grounds,
         biomass_kg_ha,
         abundance_ind_250m2) %>%
  pivot_longer(cols = c(biomass_kg_ha, abundance_ind_250m2),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         biomass_kg_ha = "Biomass (kg/ha)",
                         abundance_ind_250m2 = "Abundance (ind/250m²)"))

# Combined boxplot
ggplot(dov_long, aes(x = Fishing.Grounds, y = Value, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  facet_wrap(~ Metric, scales = "free_y") +   # separate panels for each metric
  labs(x = "Fishing Grounds", y = NULL,
       title = "Fish Biomass and Abundance by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
