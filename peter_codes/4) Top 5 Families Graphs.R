# Top 5 fish families per site by location
# Version date - 4/11/25

# Combine EmobsCombined_full and opcode_groups without duplication
Emobs_full <- EmobsCombined_Filled %>%
  left_join(opcode_groups %>% distinct(OpCode, .keep_all = TRUE),
            by = "OpCode")
head(Emobs_full)

# Compute abundance density (individuals / 250m²)
Emobs_full <- Emobs_full %>%
  mutate(abundance_ind_250m2 = 1) %>%  # each fish = 1 individual per transect
  group_by(Location, Site, Family, Period) %>%
  summarise(
    abundance_ind_250m2 = sum(abundance_ind_250m2, na.rm = TRUE),
    biomass_kg = sum(biomass_kg, na.rm = TRUE),
    .groups = "drop"
  )

# Convert to biomass density (kg/ha)
Emobs_full <- Emobs_full %>%
  mutate(biomass_kg_ha = biomass_kg * (10000 / 250))  # scale 250 m² → 1 ha

# Top 5 summaries
top_families_abund <- Emobs_full %>%
  group_by(Location, Site, Family) %>%
  summarise(mean_abundance = mean(abundance_ind_250m2, na.rm = TRUE), .groups = "drop") %>%
  group_by(Location, Site) %>%
  mutate(rank = rank(-mean_abundance, ties.method = "first")) %>%
  filter(rank <= 5)

# Biomass
top_families_biomass_ha <- Emobs_full %>%
  group_by(Location, Site, Family) %>%
  summarise(mean_biomass_ha = mean(biomass_kg_ha, na.rm = TRUE), .groups = "drop") %>%
  group_by(Location, Site) %>%
  mutate(rank = rank(-mean_biomass_ha, ties.method = "first")) %>%
  filter(rank <= 5)

#--- 3. Biomass plot ---
p1 <- ggplot(top_families_biomass_ha, aes(x = Site, y = mean_biomass_ha, fill = Family)) +
  geom_col(position = "stack", color = "black") +
  facet_wrap(~ Location, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Top 5 Fish Families by Biomass Density",
    x = "Site",
    y = "Mean Biomass (kg/ha)"
  )

#--- 4. Abundance plot ---
p2 <- ggplot(top_families_abund, aes(x = Site, y = mean_abundance, fill = Family)) +
  geom_col(position = "stack", color = "black") +
  facet_wrap(~ Location, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Top 5 Fish Families by Abundance",
    x = "Site",
    y = "Mean Abundance (individuals / 250 m²)",
    fill = "Family"
  )

#--- 5. Combine both plots side by side ---
p1 / p2 + plot_layout(heights = c(1, 1.1))

# Combine into clean figure
top5_plot <- p1 / p2 +
  plot_annotation(title = "Top 5 Families across Locations",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))
top5_plot

# Save
ggsave("Top5_by_Location.png", top5_plot,
       width = 14, height = 20, dpi = 300)

# -------------------------------------------------------------
# Top 5 fish families per site (by location)
# Version date - 5/11/25
# -------------------------------------------------------------

# Combine EmobsCombined_Filled with opcode_groups (unique by OpCode)
Emobs_full <- EmobsCombined_Filled %>%
  left_join(opcode_groups %>% distinct(OpCode, .keep_all = TRUE),
            by = "OpCode")

# Compute abundance (individuals / 250m²) and biomass
Emobs_full <- Emobs_full %>%
  mutate(abundance_ind_250m2 = 1) %>%  # each observation = 1 individual
  group_by(Location, Family, Period) %>%
  summarise(
    abundance_ind_250m2 = sum(abundance_ind_250m2, na.rm = TRUE),
    biomass_kg = sum(biomass_kg, na.rm = TRUE),
    .groups = "drop"
  )

# Convert to biomass density (kg/ha)
Emobs_full <- Emobs_full %>%
  mutate(biomass_kg_ha = biomass_kg * (10000 / 250))  # scale 250 m² → 1 ha

# -------------------------------------------------------------
# Top 5 summaries per location
# -------------------------------------------------------------

# Abundance
top_families_abund_loc <- Emobs_full %>%
  group_by(Location, Family) %>%
  summarise(mean_abundance = mean(abundance_ind_250m2, na.rm = TRUE), .groups = "drop") %>%
  group_by(Location) %>%
  mutate(rank = rank(-mean_abundance, ties.method = "first")) %>%
  filter(rank <= 5)

# Biomass
top_families_biomass_loc <- Emobs_full %>%
  group_by(Location, Family) %>%
  summarise(mean_biomass_ha = mean(biomass_kg_ha, na.rm = TRUE), .groups = "drop") %>%
  group_by(Location) %>%
  mutate(rank = rank(-mean_biomass_ha, ties.method = "first")) %>%
  filter(rank <= 5)

# -------------------------------------------------------------
# Plotting
# -------------------------------------------------------------
# --- Biomass plot ---
p1 <- ggplot(top_families_biomass_loc, aes(x = reorder(Family, -mean_biomass_ha), 
                                           y = mean_biomass_ha, fill = Family)) +
  geom_col(color = "black") +
  facet_wrap(~ Location, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Top 5 Fish Families by Biomass Density (kg/ha)",
    x = "Family",
    y = "Mean Biomass (kg/ha)"
  )

# --- Abundance plot ---
p2 <- ggplot(top_families_abund_loc, aes(x = reorder(Family, -mean_abundance),
                                         y = mean_abundance, fill = Family)) +
  geom_col(color = "black") +
  facet_wrap(~ Location, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Top 5 Fish Families by Abundance (individuals / 250 m²)",
    x = "Family",
    y = "Mean Abundance",
    fill = "Family"
  )

# --- Combine plots vertically ---
top5_plot <- p1 / p2 +
  plot_annotation(
    title = "Top 5 Fish Families across Locations",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

top5_plot

# Save
ggsave("Top5_per_Location.png", top5_plot,
       width = 14, height = 15, dpi = 300)