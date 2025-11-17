###############################################################################
# FUNCTIONAL TRAIT SPACES (Thesis Version – Fixed + Region Colours)
###############################################################################
library(ggrepel)

#-------------------------------------------------------------
# 1. Species coordinates in functional trait space (Gower + PCoA)
#-------------------------------------------------------------
gower_dist <- daisy(traits_mat2, metric = "gower")
pcoa_res <- cmdscale(gower_dist, k = 4, eig = TRUE)

point <- data.frame(
  Species_full = rownames(pcoa_res$points),
  A1 = pcoa_res$points[, 1],
  A2 = pcoa_res$points[, 2],
  A3 = pcoa_res$points[, 3],
  A4 = pcoa_res$points[, 4]
)

var_explained <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)
var_explained[1:4]

# Quick global trait check
ggplot(point, aes(x = A1, y = A2)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  labs(x = paste0("PCoA Axis 1 (", var_explained[1], "%)"),
       y = paste0("PCoA Axis 2 (", var_explained[2], "%)"),
       title = "Global Functional Trait Space of Reef Fish") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

###############################################################################
# 2. Build figure datasets (using OpCode + lookup)
###############################################################################

opcode_groups <- opcode_groups %>%
  rename(Fishing_Grounds = `Fishing Grounds`)

# Region grouping for consistent colours
opcode_groups <- opcode_groups %>%
  mutate(Region = case_when(
    grepl("Kadavu", Fishing_Grounds, ignore.case = TRUE) ~ "Kadavu",
    grepl("Namara", Fishing_Grounds, ignore.case = TRUE) ~ "Yasawas (Namara)",
    grepl("Muaira", Fishing_Grounds, ignore.case = TRUE) ~ "Yasawas (Muaira)",
    grepl("Totoya|Vanua Balavu", Fishing_Grounds, ignore.case = TRUE) ~ "Lau",
    grepl("Savusavu|Wailevu", Fishing_Grounds, ignore.case = TRUE) ~ "Vanua Levu",
    TRUE ~ "Other"
  ))

# 2.1 OpCode × Species abundance
fish_assemblage <- EmobsCombined_Filled %>%
  group_by(OpCode, Species_full) %>%
  summarise(abund = n(), .groups = "drop") %>%
  left_join(opcode_groups %>% select(OpCode, Site, Fishing_Grounds, Region),
            by = "OpCode")

# 2.2 Total abundance per OpCode
site.total <- fish_assemblage %>%
  group_by(OpCode) %>%
  summarise(site_total = sum(abund), .groups = "drop") %>%
  left_join(opcode_groups %>% select(OpCode, Site, Fishing_Grounds, Region),
            by = "OpCode")

# 2.3 Trait table for plotting
traits <- traits_mat2 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species_full") %>%
  rename(ForageMode = Diet_Mouillot_2014)

# 2.4 FD metrics
transect_meta <- EmobsCombined_Filled %>%
  distinct(Period, OpCode)

fd_site <- fd_metrics %>%
  left_join(transect_meta, by = c("Transect" = "Period")) %>%
  group_by(OpCode) %>%
  summarise(
    TOP  = mean(FRic, na.rm = TRUE),
    TEve = mean(FEve, na.rm = TRUE),
    TDiv = mean(FDiv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(opcode_groups %>% select(OpCode, Site, Fishing_Grounds, Region),
            by = "OpCode")

# 2.5 Functional group (diet) composition
func_groups <- EmobsCombined_Filled %>%
  left_join(combined_traits %>% select(Species_full, Diet_Mouillot_2014),
            by = "Species_full") %>%
  left_join(opcode_groups %>% select(OpCode, Site, Fishing_Grounds, Region),
            by = "OpCode") %>%
  group_by(Site, Diet_Mouillot_2014) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(rel_abund = n / sum(n)) %>%
  pivot_wider(id_cols = Site,
              names_from = Diet_Mouillot_2014,
              values_from = rel_abund,
              values_fill = 0)

predictors <- left_join(fd_site, func_groups, by = "Site")

###############################################################################
# 3. Helper functions and colours
###############################################################################

convhull.vert <- function(df) df[chull(df), , drop = FALSE]

# Region colour palette
region_cols <- c(
  "Kadavu" = "#7B68EE",
  "Yasawas (Namara)" = "#D62728",
  "Yasawas (Muaira)" = "#FFD700",
  "Lau" = "#FF8C00",
  "Vanua Levu" = "#1F77B4"
)

###############################################################################
# 4. Build per-site convex hull plots
###############################################################################

titles <- unique(fish_assemblage$Site)
t.space <- vector("list", length(titles))

hull.v  <- convhull.vert(point[, 2:3])
hull.v2 <- convhull.vert(point[, 4:5])

for (i in seq_along(titles)) {
  site_name <- titles[i]
  
  site_row <- site.total %>% filter(Site == site_name)
  if (nrow(site_row) == 0) next
  
  site_fg <- unique(site_row$Fishing_Grounds)
  site_region <- unique(site_row$Region)
  site_col <- region_cols[site_region]
  
  pointS <- fish_assemblage %>%
    filter(Site == site_name) %>%
    left_join(point, by = "Species_full") %>%
    mutate(rel_abund = abund / site_row$site_total)
  
  if (nrow(pointS) == 0) next
  
  s.hull.v  <- convhull.vert(pointS[, c("A1","A2")])
  s.hull.v2 <- convhull.vert(pointS[, c("A3","A4")])
  
  top_species <- pointS %>% arrange(desc(rel_abund)) %>% slice_head(n = 3)
  shortsp <- sapply(strsplit(top_species$Species_full, " "),
                    function(x) paste0(substr(x[1], 1, 1), ". ", x[2]))
  pointS$ShortSp <- ""
  pointS$ShortSp[match(top_species$Species_full, pointS$Species_full)] <- shortsp
  
  p1 <- ggplot() + theme_bw(base_size = 11) +
    geom_polygon(data = hull.v, aes(x = A1, y = A2), color = 'grey70', fill = 'white') +
    geom_polygon(data = s.hull.v, aes(x = A1, y = A2), fill = site_col, alpha = 0.4) +
    geom_point(data = pointS, aes(x = A1, y = A2, size = rel_abund),
               shape = 21, fill = site_col, color = "grey30") +
    geom_text_repel(data = top_species, aes(x = A1, y = A2, label = shortsp),
                    size = 2.3, fontface = "italic", box.padding = 0.5) +
    labs(x = "PCo1", y = "PCo2",
         title = paste0(site_name, "\n(", site_fg, ")")) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_size_continuous(range = c(1.5, 8), guide = "none")
  
  p2 <- ggplot() + theme_bw(base_size = 11) +
    geom_polygon(data = hull.v2, aes(x = A3, y = A4), color = 'grey70', fill = 'white') +
    geom_polygon(data = s.hull.v2, aes(x = A3, y = A4), fill = site_col, alpha = 0.4) +
    geom_point(data = pointS, aes(x = A3, y = A4, size = rel_abund),
               shape = 21, fill = site_col, color = "grey30") +
    geom_text_repel(data = top_species, aes(x = A3, y = A4, label = shortsp),
                    size = 2.3, fontface = "italic", box.padding = 0.5) +
    labs(x = "PCo3", y = "PCo4") +
    theme(panel.grid = element_blank()) +
    scale_size_continuous(range = c(1.5, 8), guide = "none")
  
  t.space[[i]] <- p1 / p2
}

###############################################################################
# 5. Global trait space (all sites combined)
###############################################################################

Sp.count <- fish_assemblage %>%
  group_by(Species_full) %>%
  summarise(n = sum(abund), .groups="drop")

point$n <- Sp.count$n[match(point$Species_full, Sp.count$Species_full)]
point$n <- point$n / sum(point$n, na.rm = TRUE)

global1 <- ggplot() +
  geom_polygon(data = hull.v, aes(x = A1, y = A2), fill = "transparent", color = "grey") +
  geom_point(data = left_join(point, traits, by = "Species_full"),
             aes(x = A1, y = A2, shape = ForageMode), size = 2, color = "black", alpha = 0.6) +
  labs(x = "PCo1", y = "PCo2", title = "All sites") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_shape_manual(values = 1:10, name = "Foraging mode")

global2 <- ggplot() +
  geom_polygon(data = hull.v2, aes(x = A3, y = A4), fill = "transparent", color = "grey") +
  geom_point(data = left_join(point, traits, by = "Species_full"),
             aes(x = A3, y = A4, shape = ForageMode), size = 2, color = "black", alpha = 0.6) +
  labs(x = "PCo3", y = "PCo4") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank()) +
  scale_shape_manual(values = 1:10, guide = "none")

###############################################################################
# 6. Combine final figure
###############################################################################

bars <- (ggplot() + theme_void())  # placeholder if you already built barplots earlier

sites_panel <- wrap_plots(t.space, ncol = 4)

Fig2 <- (global1 / global2 / sites_panel) *
  plot_layout(heights = c(1.2, 1.2, 8)) &
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

###############################################################################
# 7. Export (A3 thesis quality)
###############################################################################

ggsave("TraitSpaces_Thesis_Final.svg",
       plot = Fig2,
       width = 42, height = 60, units = "cm", dpi = 600)
###############################################################################