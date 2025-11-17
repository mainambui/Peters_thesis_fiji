# Biomass and Abundance Graphs
# Site level
# Version date - 4/11/25

# Mean biomass and abundance
Sum_Biomass <- read.csv("Sum_Biomass.csv")
Sum_Abundance <- read.csv("Sum_Abundance.csv")
opcode_groups <- read.csv("opcode_groups.csv")

# Combine the two summaries by OpCode (transect)
combined_df <- Sum_Biomass %>%
  left_join(Sum_Abundance, by = c("OpCode", "Period")) %>%
  left_join(opcode_groups, by = "OpCode") %>%
  select(
    Period, OpCode,
    biomass_kg, biomass_kg_ha,
    total_abundance, abundance_ind_250m2,
    Location, Site, Year.of.Adoption, Status, Province
  )

# Summarize by site and location
summary_df <- combined_df %>%
  group_by(Location, Site) %>%
  summarise(
    mean_biomass = mean(biomass_kg_ha, na.rm = TRUE),
    se_biomass = sd(biomass_kg_ha, na.rm = TRUE) / sqrt(sum(!is.na(biomass_kg_ha))),
    .groups = "drop"
  )

# Biomass
ggplot(summary_df, aes(x = Site, y = mean_biomass, fill = Location)) +
  geom_col(position = "dodge") +
  geom_errorbar(
    aes(ymin = mean_biomass - se_biomass, ymax = mean_biomass + se_biomass),
    width = 0.2
  ) +
  facet_wrap(~ Location, scales = "free_x", ncol = 2) +  # 4 locations, 2 columns layout
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",  # color already shown by facet
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = "Mean Fish Biomass by Site and Location",
    x = "Site",
    y = "Biomass (kg/ha)"
  )

# Abundance
summary_abund <- combined_df %>%
  group_by(Location, Site) %>%
  summarise(
    mean_abundance = mean(abundance_ind_250m2, na.rm = TRUE),
    se_abundance = sd(abundance_ind_250m2, na.rm = TRUE) / sqrt(sum(!is.na(abundance_ind_250m2))),
    .groups = "drop"
  )

ggplot(summary_abund, aes(x = Site, y = mean_abundance, fill = Location)) +
  geom_col(position = "dodge") +
  geom_errorbar(
    aes(ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance),
    width = 0.2
  ) +
  facet_wrap(~ Location, scales = "free_x", ncol = 2) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = "Mean Fish Abundance by Site and Location",
    x = "Site",
    y = "Abundance (individuals / 250 m²)"
  )

# Sites ordered consistently between both plots
# Name plots
# ---- Biomass plot ----
p_bio <- ggplot(summary_df, aes(x = Site, y = mean_biomass, fill = Location)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = mean_biomass - se_biomass, ymax = mean_biomass + se_biomass),
                width = 0.2) +
  facet_wrap(~ Location, scales = "free_x", ncol = 2) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = "Mean Fish Biomass by Site and Location",
    x = "Site",
    y = "Biomass (kg/ha)"
  )

# ---- Abundance plot ----
p_abund <- ggplot(summary_abund, aes(x = Site, y = mean_abundance, fill = Location)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance),
                width = 0.2) +
  facet_wrap(~ Location, scales = "free_x", ncol = 2) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = "Mean Fish Abundance by Site and Location",
    x = "Site",
    y = "Abundance (individuals / 250 m²)"
  )

# Combine into clean figure
combined_plot <- p_bio / p_abund +
  plot_annotation(title = "Fish Biomass and Abundance across Locations",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))
combined_plot

# Save
ggsave("Biomass_Abundance_by_Location.png", combined_plot,
       width = 14, height = 10, dpi = 300)

# Biomass and Abundance Graphs
# Version date - 5/11/25
# Load data
Sum_Biomass <- read.csv("Sum_Biomass.csv") %>% mutate(OpCode = as.character(OpCode))
Sum_Abundance <- read.csv("Sum_Abundance.csv") %>% mutate(OpCode = as.character(OpCode))
opcode_groups <- read.csv("opcode_groups.csv") %>% mutate(OpCode = as.character(OpCode))

# Combine and summarise per site
summary_combined <- Sum_Biomass %>%
  left_join(Sum_Abundance, by = c("OpCode", "Period")) %>%
  left_join(opcode_groups, by = "OpCode") %>%
  group_by(Location, Site) %>%
  summarise(
    mean_biomass = mean(biomass_kg_ha, na.rm = TRUE),
    se_biomass = sd(biomass_kg_ha, na.rm = TRUE) / sqrt(sum(!is.na(biomass_kg_ha))),
    mean_abundance = mean(abundance_ind_250m2, na.rm = TRUE),
    se_abundance = sd(abundance_ind_250m2, na.rm = TRUE) / sqrt(sum(!is.na(abundance_ind_250m2))),
    .groups = "drop"
  ) %>%
  arrange(Location, Site) %>%
  mutate(Site = factor(Site, levels = unique(Site)))

# --- Pivot to long format for side-by-side plotting ---
summary_long <- summary_combined %>%
  pivot_longer(
    cols = c(mean_biomass, mean_abundance),
    names_to = "Metric",
    values_to = "Mean"
  ) %>%
  mutate(
    SE = case_when(
      Metric == "mean_biomass" ~ se_biomass,
      Metric == "mean_abundance" ~ se_abundance
    ),
    Metric = recode(Metric,
                    "mean_biomass" = "Biomass (kg/ha)",
                    "mean_abundance" = "Abundance (ind / 250 m²)")
  )

# --- Plot ---
ggplot(summary_long, aes(x = Site, y = Mean, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  facet_wrap(~ Location, scales = "free_x", ncol = 2) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    title = "Fish Biomass and Abundance by Site and Location",
    x = "Site",
    y = "Value (kg/ha or individuals / 250 m²)",
    caption = "Error bars = SE"
  )

# --- Save ---
if (!dir.exists("figures")) dir.create("figures")
ggsave("figures/Biomass_Abundance_SideBySide_ByLocation.png",
       width = 14, height = 8, dpi = 400)