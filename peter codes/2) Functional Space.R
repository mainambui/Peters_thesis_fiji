# Run after transect  response variables
names(fd_results)
#------------------------------------------------------------
#--Extract the PCoA axes manually (independent of FD)--
#------------------------------------------------------------
# Step 1: Compute Gower distance matrix from your trait matrix
library(cluster)
gower_dist <- daisy(traits_mat2, metric = "gower")

# Step 2: Run PCoA using Cailliez correction (same as dbFD)
library(vegan)
pcoa_res <- cmdscale(as.dist(gower_dist), k = 4, add = TRUE)  

# Step 3: Extract coordinates
species_axes <- as.data.frame(pcoa_res$points)

# Add species names
species_axes <- species_axes %>%
  rownames_to_column("Species_full")

# Rename PCoA axes
colnames(species_axes)[colnames(species_axes) == "V1"] <- "A1"
colnames(species_axes)[colnames(species_axes) == "V2"] <- "A2"
colnames(species_axes)[colnames(species_axes) == "V3"] <- "A3"
colnames(species_axes)[colnames(species_axes) == "V4"] <- "A4"

# Clean final order
species_axes <- species_axes %>%
  select(Species_full, A1, A2, A3, A4)

# Sep 4: Merge with traits if needed
species_axes_traits <- species_axes %>%
  left_join(combined_traits, by = "Species_full")

#--------------------------------------------------------------
#global + per–Fishing Grounds trait spaces using:
#PCoA axes (species_axes), species traits, and abundance data.
#--------------------------------------------------------------
# STEP 1 — Build abundance per Fishing Ground × Species
fish_fg <- EmobsCombined_Filled %>%
  left_join(opcode_groups, by = "OpCode") %>%   # add location, fishing grounds, status, etc.
  group_by(`Fishing Grounds`, Species_full) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(`Fishing Grounds`) %>%
  mutate(rel_abund = n / sum(n)) %>%
  ungroup()
#Check fishing grounds present
unique(fish_fg$`Fishing Grounds`)

#------------------------------------------------------------
# REMOVE species with missing PCoA axes (cannot plot hulls)
#------------------------------------------------------------
species_axes_clean <- species_axes_traits %>%
  filter(!is.na(A1), !is.na(A2), !is.na(A3), !is.na(A4))

# Clean fish_fg to include only species with PCoA coords
fish_fg <- fish_fg %>%
  inner_join(species_axes_clean %>% select(Species_full), by = "Species_full")

# Remove one fishing ground
fish_fg <- fish_fg %>% 
  filter(`Fishing Grounds` != "Vanua Balavu")

# STEP 2 — Build global convex hull (A1–A2, A3–A4)
# convex hull function
convhull.vert <- function(df) {
  df <- as.data.frame(df)
  df[] <- lapply(df, as.numeric)
  
  df <- df[complete.cases(df), , drop = FALSE]
  
  if (nrow(df) < 3) {
    return(NULL)   # THIS solves the error
  }
  
  h <- chull(df[,1], df[,2])
  df[h, , drop = FALSE]
}

# GLOBAL hulls
hull12 <- convhull.vert(species_axes_clean[, c("A1","A2")])
hull34 <- convhull.vert(species_axes_clean[, c("A3","A4")])

# STEP 3 — Create GLOBAL TRAIT SPACE plot
library(ggplot2)
library(ggrepel)
library(viridis)

theme_trait <- theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.3))

globalA <- ggplot() +
  theme_trait +
  geom_polygon(data = hull12, aes(A1, A2),
               fill = NA, color = "black", linewidth = 1) +
  geom_point(data = species_axes_clean,
             aes(A1, A2, shape = Diet_Mouillot_2014),
             size = 2, fill = "black", alpha = 0.7) +
  labs(title = "Global Trait Space: PCoA1 × PCoA2",
       x = "PCoA1", y = "PCoA2") +
  scale_shape_manual(values = 1:10)

globalB <- ggplot() +
  theme_trait +
  geom_polygon(data = hull34, aes(A3, A4),
               fill = NA, color = "black", linewidth = 1) +
  geom_point(data = species_axes_clean,
             aes(A3, A4, shape = Diet_Mouillot_2014),
             size = 2, fill = "black", alpha = 0.7) +
  labs(title = "Global Trait Space: PCoA3 × PCoA4",
       x = "PCoA3", y = "PCoA4") +
  scale_shape_manual(values = 1:10)

# STEP 4 — Plot per–Fishing Ground trait spaces
fg_list <- unique(fish_fg$`Fishing Grounds`)
n_fg <- length(fg_list)

common_size_scale <- scale_size_continuous(
  range = c(2, 8),
  name = "Relative abundance")

fg_plots <- vector("list", n_fg)
palette_fg <- viridis(n_fg)

for (i in seq_len(n_fg)) {
  
  fg_i <- fg_list[i]
  
  df_i <- fish_fg %>%
    filter(`Fishing Grounds` == fg_i) %>%
    left_join(species_axes_clean, by = "Species_full") %>%
    arrange(desc(rel_abund)) %>%
    mutate(
      ShortSp = case_when(
        row_number() <= 3 ~ paste0(substr(word(Species_full,1),1,1), ". ", word(Species_full,2)),
        TRUE ~ ""
      )
    )
  
  # hulls
  hull_i_12 <- convhull.vert(df_i[, c("A1","A2")])
  hull_i_34 <- convhull.vert(df_i[, c("A3","A4")])
  
  # TOP panel (keeps legend)
  p12 <- ggplot() +
    theme_trait +
    geom_polygon(data = hull12, aes(A1, A2), fill="white", color="grey60") +
    geom_polygon(data = hull_i_12, aes(A1, A2),
                 fill = palette_fg[i], alpha = 0.4) +
    geom_point(data = df_i,
               aes(A1, A2, size = rel_abund),
               shape = 21, fill = palette_fg[i], color = "black") +
    geom_text_repel(data = df_i %>% filter(ShortSp != ""),
                    aes(A1, A2, label = ShortSp), size = 3) +
    labs(title = fg_i, x = "PCoA1", y = "PCoA2") +
    common_size_scale +
    guides(size = ifelse(i == 1, "legend", "none"))   # <-- THIS IS THE FIX
  
  # BOTTOM panel (no legend)
  p34 <- ggplot() +
    theme_trait +
    geom_polygon(data = hull34, aes(A3, A4), fill="white", color="grey60") +
    geom_polygon(data = hull_i_34, aes(A3, A4),
                 fill = palette_fg[i], alpha = 0.4) +
    geom_point(data = df_i,
               aes(A3, A4, size = rel_abund),
               shape = 21, fill = palette_fg[i], color = "black") +
    geom_text_repel(data = df_i %>% filter(ShortSp != ""),
                    aes(A3, A4, label = ShortSp), size = 3) +
    labs(x = "PCoA3", y = "PCoA4") +
    common_size_scale +
    guides(size = "none")     # <--- THIS IS THE CRITICAL FIX
  
  # Combine into a 2-row FG plot
  fg_plots[[i]] <- p12 / p34
}

#------------------------------------------------------------------------------
# FINAL FIGURE CODE
#------------------------------------------------------------------------------
library(patchwork)

# Ensure axis limits are consistent across all FG plots
x1_lim <- range(species_axes_clean$A1, na.rm = TRUE)
x2_lim <- range(species_axes_clean$A2, na.rm = TRUE)
x3_lim <- range(species_axes_clean$A3, na.rm = TRUE)
x4_lim <- range(species_axes_clean$A4, na.rm = TRUE)

# Add fixed limits to global plots
globalA <- globalA +
  coord_cartesian(xlim = x1_lim, ylim = x2_lim)

globalB <- globalB +
  coord_cartesian(xlim = x3_lim, ylim = x4_lim)

# Add fixed limits to all FG plots
fg_plots_fixed <- lapply(fg_plots, function(p){
  p + 
    coord_cartesian(expand = FALSE) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
})

#-------------------------------------------
#  FIGURE 5 — Style 2: Condensed Journal Layout
#-------------------------------------------

library(patchwork)

# --- Create a blank legend placeholder ---
legend_placeholder <- ggplot() + 
  theme_void() +
  ggtitle("Legend") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# --- Top row: Global A1–A2 | Global A3–A4 | Legend ---
top_row <- globalA | globalB | legend_placeholder
top_row <- top_row + plot_layout(widths = c(1, 1, 0.6), guides = "collect")

# --- 3-column FG grid ---
fg_grid <- wrap_plots(fg_plots_fixed, ncol = 3) +
  plot_layout(guides = "collect")

# --- Final Figure 5 assembly ---
TraitSpace_All <- top_row /
  fg_grid +
  plot_layout(
    heights = c(1, 3),
    guides = "collect"
  ) +
  plot_annotation(
    title = "Figure 5. Functional Trait Space — Global and by Fishing Ground",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.margin = margin(6, 6, 6, 6)
    )
  )

# Save
ggsave("TraitSpace_FishingGrounds.png",
       plot = TraitSpace_All,
       width = 14,
       height = 16,
       dpi = 300,
       bg = "white")