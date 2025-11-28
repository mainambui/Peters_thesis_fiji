# Run after 4) Top 5 families
# ================================================================
#  Load data
# ================================================================
dov_final <- read.csv("DOV_Transects_Properties_With_Biomass_Abundance_FD.csv")

# ================================================================
#  Exclude Daliconi & Prepare Base Dataset
# ================================================================
dov_plot <- dov_final %>%
  filter(Fishing.Grounds != "Daliconi") %>%
  select(Fishing.Grounds,
         biomass_kg_ha,
         abundance_ind_250m2,
         FRic, FEve, FDiv, FDis, RaoQ)

# ================================================================
#  1. Combined Biomass + Abundance Boxplot (long format)
# ================================================================
dov_long <- dov_final %>%
  filter(Fishing.Grounds != "Daliconi") %>%
  select(Fishing.Grounds, biomass_kg_ha, abundance_ind_250m2) %>%
  pivot_longer(cols = c(biomass_kg_ha, abundance_ind_250m2),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         biomass_kg_ha = "Biomass (kg/ha)",
                         abundance_ind_250m2 = "Abundance (ind/250m²)"))

combined_biomass_abundance <- ggplot(dov_long,
                                     aes(x = Fishing.Grounds, y = Value, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(x = "Fishing Grounds", y = NULL,
       title = "Fish Biomass and Abundance by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display
combined_biomass_abundance

# Save
#ggsave("Combined_Biomass_Abundance.png",
#       combined_biomass_abundance, dpi = 300, width = 8, height = 5)


# ================================================================
# 2. Biomass Plots (3 versions)
# ================================================================

# --- Biomass: Boxplot ---
biomass_box <- ggplot(dov_plot,
                      aes(x = Fishing.Grounds, y = biomass_kg_ha, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  labs(x = "Fishing Grounds", y = "Fish Biomass (kg/ha)",
       title = "Fish Biomass by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display
biomass_box

# Save
#ggsave("Biomass_Boxplot.png",
#       biomass_box, dpi = 300, width = 7, height = 5)


# --- Biomass: Boxplot + Jitter + Mean ---
biomass_jitter_mean <- ggplot(dov_plot,
                              aes(x = Fishing.Grounds, y = biomass_kg_ha, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 23,
               size = 4, fill = "red", color = "black") +
  labs(x = "Fishing Grounds", y = "Fish Biomass (kg/ha)",
       title = "Biomass") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display
biomass_jitter_mean

# Save
#ggsave("Biomass_Boxplot_Jitter_Mean.png",
#       biomass_jitter_mean, dpi = 300, width = 7, height = 5)


# --- Biomass: Boxplot + Jitter + Mean ± SE ---
biomass_mean_se <- ggplot(dov_plot,
                          aes(x = Fishing.Grounds, y = biomass_kg_ha, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 23,
               size = 4, fill = "red", color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.3, color = "red", size = 0.8) +
  labs(x = "Fishing Grounds", y = "Fish Biomass (kg/ha)",
       title = "Fish Biomass by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

#Display
biomass_mean_se

# Save
#ggsave("Biomass_Mean_SE.png",
#       biomass_mean_se, dpi = 300, width = 7, height = 5)


# ================================================================
# 3. Abundance Plots (3 versions)
# ================================================================

# --- Abundance: Boxplot ---
abundance_box <- ggplot(dov_plot,
                        aes(x = Fishing.Grounds, y = abundance_ind_250m2, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  labs(x = "Fishing Grounds", y = "Fish Abundance (individuals / 250m²)",
       title = "Fish Abundance by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display
abundance_box

# Save
#ggsave("Abundance_Boxplot.png",
#       abundance_box, dpi = 300, width = 7, height = 5)


# --- Abundance: Boxplot + Jitter + Mean ---
abundance_jitter_mean <- ggplot(dov_plot,
                                aes(x = Fishing.Grounds, y = abundance_ind_250m2, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 23,
               size = 4, fill = "red", color = "black") +
  labs(x = "Fishing Grounds",
       y = "Fish Abundance (individuals / 250m²)",
       title = "Abundance") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display
abundance_jitter_mean

# Save
#ggsave("Abundance_Boxplot_Jitter_Mean.png",
#       abundance_jitter_mean, dpi = 300, width = 7, height = 5)


# --- Abundance: Boxplot + Jitter + Mean ± SE ---
abundance_mean_se <- ggplot(dov_plot,
                            aes(x = Fishing.Grounds, y = abundance_ind_250m2, fill = Fishing.Grounds)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 23,
               size = 4, fill = "red", color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.3, color = "red", size = 0.8) +
  labs(x = "Fishing Grounds",
       y = "Fish Abundance (individuals / 250m²)",
       title = "Fish Abundance by Fishing Grounds") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display
abundance_mean_se

#Save
#ggsave("Abundance_Mean_SE.png",
#       abundance_mean_se, dpi = 300, width = 7, height = 5)


# ================================================================
# 4. Combined Biomass × Abundance Patchwork
# ================================================================
#combined_ba <- (biomass_jitter_mean | abundance_jitter_mean) +
#  plot_annotation(title = "Fish Biomass and Abundance by Fishing Grounds",
#                  theme = theme(plot.title = element_text(size = 18, face = "bold")))
# Display
#combined_ba
#ggsave("Biomass_Abundance_Combined.png",
#       combined_ba, dpi = 300, width = 8, height = 10)

# ================================================================
# Combined Biomass × Abundance (facet_wrap version)
# ================================================================

combined_ba_facet <- ggplot(dov_long,
                            aes(x = Fishing.Grounds, y = Value, fill = Fishing.Grounds)) +
  # boxplot (semi-transparent)
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  
  # jitter points
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, colour = "black") +
  
  # mean points
  stat_summary(
    fun = mean, geom = "point", shape = 23, size = 4,
    fill = "red", colour = "black"
  ) +
  
  # facet layout
  facet_wrap(~ Metric, scales = "free_y") +
  
  # labels
  labs(
    x = "Fishing Grounds",
    y = NULL,
    title = "Fish Biomass and Abundance by Fishing Grounds"
  ) +
  
  # theme
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# Display
combined_ba_facet

# Save
ggsave("Biomass_Abundance_Combined.png",
       combined_ba_facet, dpi = 300, width = 8, height = 6)

# ================================================================
# 5. Functional Diversity (FD) Individual Metrics WITH transect points + mean ONLY
# ================================================================

make_fd_plot <- function(metric, label) {
  ggplot(dov_plot,
         aes_string(x = "Fishing.Grounds", y = metric, fill = "Fishing.Grounds")) +
    
    # Boxplot
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    
    # Transect-level points
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +
    
    # Mean per fishing ground (no SE bars)
    stat_summary(fun = mean, geom = "point",
                 shape = 23, size = 4,
                 fill = "red", color = "black") +
    
    # Labels & theme
    labs(x = "Fishing Grounds", y = label,
         title = paste(label, "by Fishing Grounds")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

fd_metrics <- list(
  FRic = "Functional Richness (FRic)",
  FEve = "Functional Evenness (FEve)",
  FDiv = "Functional Divergence (FDiv)",
  FDis = "Functional Dispersion (FDis)",
  RaoQ = "Rao's Quadratic Entropy (RaoQ)"
)

for (m in names(fd_metrics)) {
  p <- make_fd_plot(m, fd_metrics[[m]])
  print(p)  # optional preview
  fname <- paste0("FD_", m, ".png")
  ggsave(fname, p, dpi = 300, width = 7, height = 5)
}

# ================================================================
# 6. Combined FD Facet Plot (NO SE bars)
# ================================================================

dov_fd_long <- dov_plot %>%
  pivot_longer(cols = c(FRic, FEve, FDiv, FDis, RaoQ),
               names_to = "FD_metric",
               values_to = "value")

fd_facet <- ggplot(dov_fd_long,
                   aes(x = Fishing.Grounds, y = value, fill = Fishing.Grounds)) +
  
  # Boxplot
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  
  # Transect-level points
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "black") +
  
  # Mean per fishing ground (no SE bars)
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 3,
               fill = "red", color = "black") +
  
  facet_wrap(~ FD_metric, scales = "free_y") +
  
  labs(x = "Fishing Grounds",
       y = "Functional Diversity Metric",
       title = "Functional Diversity Metrics by Fishing Grounds") +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display
fd_facet

# Save
ggsave("FD_Facet_All_Metrics.png",
       fd_facet, dpi = 300, width = 10, height = 8)