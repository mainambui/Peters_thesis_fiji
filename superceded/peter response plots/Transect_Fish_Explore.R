# Load packages
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(patchwork)

# Read file 
dov_final <- read.csv("DOV_transect&connect_final.csv")
glimpse(dov_final)

# Prepare data
dov_plot <- dov_final %>%
  select(Fishing.Grounds,
         biomass_kg_ha,
         abundance_ind_250m2,
         FRic, FEve, FDiv, FDis, RaoQ)

# Biomass plot
ggplot(dov_plot, aes(x = Fishing.Grounds, y = biomass_kg_ha, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  labs(x = "Fishing Grounds", y = "Fish Biomass (kg/ha)", title = "Fish Biomass by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot + jittered points + mean point
ggplot(dov_plot, aes(x = Fishing.Grounds, y = biomass_kg_ha, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  # lighter fill, hide outlier dots
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # each transect
  stat_summary(fun = mean, geom = "point", 
               shape = 23, size = 4, fill = "red", color = "black") +  # mean per group
  labs(x = "Fishing Grounds",
       y = "Fish Biomass (kg/ha)",
       title = "Fish Biomass by Fishing Grounds (with Transect Data and Mean)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot + Transect Points + Mean ± SE Bars
ggplot(dov_plot, aes(x = Fishing.Grounds, y = biomass_kg_ha, fill = Fishing.Grounds)) +
  # Boxplot for distribution
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  
  # Individual transect points (each observation)
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +
  
  # Mean point per fishing ground
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 4, fill = "red", color = "black") +
  
  # Mean ± SE error bars
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.3, color = "red", size = 0.8) +
  
  # Labels and theme
  labs(x = "Fishing Grounds",
       y = "Fish Biomass (kg/ha)",
       title = "Fish Biomass by Fishing Grounds (with Transect Data, Mean ± SE)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Abundance plot
ggplot(dov_plot, aes(x = Fishing.Grounds, y = abundance_ind_250m2, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  labs(x = "Fishing Grounds", y = "Fish Abundance (individuals / 250m²)", title = "Fish Abundance by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot + jittered points + mean point
ggplot(dov_plot, aes(x = Fishing.Grounds, y = abundance_ind_250m2, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  # lighter fill, hide outlier dots
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # each transect
  stat_summary(fun = mean, geom = "point", 
               shape = 23, size = 4, fill = "red", color = "black") +  # mean per group
  labs(x = "Fishing Grounds",
       y = "Fish Abundance (individuals / 250m²)",
       title = "Fish Abundance by Fishing Grounds (with Transect Data and Mean)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot + Transect Points + Mean ± SE Bars
ggplot(dov_plot, aes(x = Fishing.Grounds, y = abundance_ind_250m2, fill = Fishing.Grounds)) +
  # Boxplot for distribution
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  
  # Individual transect points (each observation)
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +
  
  # Mean point per fishing ground
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 4, fill = "red", color = "black") +
  
  # Mean ± SE error bars
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.3, color = "red", size = 0.8) +
  
  # Labels and theme
  labs(x = "Fishing Grounds",
       y = "Fish Abundance (individuals / 250m²)",
       title = "Fish Abundance by Fishing Grounds (with Transect Data, Mean ± SE)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Biomass x Abundance 
# Biomass plot
biomass_plot <- ggplot(dov_plot, aes(x = Fishing.Grounds, y = biomass_kg_ha, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 4, fill = "red", color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.3, color = "red", size = 0.8) +
  labs(x = "Fishing Grounds",
       y = "Fish Biomass (kg/ha)",
       title = "Fish Biomass by Fishing Grounds (with Transect Data, Mean ± SE)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Abundance plot
abundance_plot <- ggplot(dov_plot, aes(x = Fishing.Grounds, y = abundance_ind_250m2, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 4, fill = "red", color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.3, color = "red", size = 0.8) +
  labs(x = "Fishing Grounds",
       y = "Fish Abundance (individuals / 250m²)",
       title = "Fish Abundance by Fishing Grounds (with Transect Data, Mean ± SE)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Combined using patchwork, stacked:
# Vertically
biomass_plot / abundance_plot
# Side by side
#biomass_plot | abundance_plot
# Add shared title/annotation
(biomass_plot / abundance_plot) +
  plot_annotation(title = "Fish Biomass and Abundance by Fishing Grounds",
                  theme = theme(plot.title = element_text(size = 18, face = "bold")))

# Richness
ggplot(dov_plot, aes(x = Fishing.Grounds, y = FRic, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Fishing Grounds", y = "FRic", title = "Functional Richness (FRic) by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Evenness
ggplot(dov_plot, aes(x = Fishing.Grounds, y = FEve, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Fishing Grounds", y = "FEve", title = "Functional Evenness (FEve) by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Divergence
ggplot(dov_plot, aes(x = Fishing.Grounds, y = FDiv, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Fishing Grounds", y = "FDiv", title = "Functional Divergence (FDiv) by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Dispersion
ggplot(dov_plot, aes(x = Fishing.Grounds, y = FDis, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Fishing Grounds", y = "FDis", title = "Functional Dispersion (FDis) by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# RaoQ
ggplot(dov_plot, aes(x = Fishing.Grounds, y = RaoQ, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Fishing Grounds", y = "RaoQ", title = "Rao's Quadratic Entropy (RaoQ) by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Combined Facet Plot - FD
dov_fd_long <- dov_plot %>%
  pivot_longer(cols = c(FRic, FEve, FDiv, FDis, RaoQ),
               names_to = "FD_metric",
               values_to = "value")

ggplot(dov_fd_long, aes(x = Fishing.Grounds, y = value, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  facet_wrap(~ FD_metric, scales = "free_y") +
  labs(x = "Fishing Grounds", y = "Functional Diversity Metric",
       title = "Functional Diversity Metrics by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90"),
        axis.text.x = element_text(angle = 45, hjust = 1))