# Must immediately run after FD calulations
#-------------------------------------------------------------
# Create species coordinates in functional trait space manually
#-------------------------------------------------------------
# Compute a Gower distance matrix from your trait matrix
gower_dist <- daisy(traits_mat2, metric = "gower")

# Run a Principal Coordinates Analysis (PCoA)
pcoa_res <- cmdscale(gower_dist, k = 4, eig = TRUE)

# Extract first 4 axes as "A1–A4"
point <- data.frame(
  Species_full = rownames(pcoa_res$points),
  A1 = pcoa_res$points[, 1],
  A2 = pcoa_res$points[, 2],
  A3 = pcoa_res$points[, 3],
  A4 = pcoa_res$points[, 4]
)

# Check result
head(point)

# Proportion of variance explained by first 4 axes
var_explained <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)
var_explained[1:4]

# Quick check: Global trait space (axes 1 & 2)
ggplot(point, aes(x = A1, y = A2)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  labs(x = paste0("PCoA Axis 1 (", var_explained[1], "%)"),
       y = paste0("PCoA Axis 2 (", var_explained[2], "%)"),
       title = "Global Functional Trait Space of Reef Fish") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

###############################################################################
# FIGURE 2 – Functional Trait Space of Reef Fish Assemblages (Thesis Version)
# Includes renamed sites and grouped by Fishing Grounds
###############################################################################
library(ggrepel)

theme_clean <- theme_bw(base_size=11) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color="grey70", size=0.3),
        axis.ticks = element_line(size=0.2),
        plot.title = element_text(face="bold", hjust=0.5, size=10))

###############################################################################
# 1. Build figure datasets (using OpCode + lookup)
###############################################################################

# 1.1 OpCode × Species abundance
fish_assemblage <- EmobsCombined_Filled %>%
  group_by(OpCode, Species_full) %>%
  summarise(n = n(), .groups="drop") %>%
  left_join(opcode_groups %>% select(OpCode, Site, Fishing_Grounds),
            by = "OpCode")

# 1.2 Total abundance per OpCode
site.total <- fish_assemblage %>%
  group_by(OpCode) %>%
  summarise(site.total = sum(n), .groups="drop") %>%
  left_join(opcode_groups %>% select(OpCode, Site, Fishing_Grounds),
            by = "OpCode")

# 1.3 Trait table for plotting
traits <- traits_mat2 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species_full") %>%
  rename(ForageMode = Diet_Mouillot_2014)

# 1.4 FD site-level metrics + functional group abundances
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
  left_join(opcode_groups %>% select(OpCode, Site, Fishing_Grounds),
            by = "OpCode")

func_groups <- EmobsCombined_Filled %>%
  left_join(combined_traits %>% select(Species_full, Diet_Mouillot_2014),
            by = "Species_full") %>%
  left_join(opcode_groups %>% select(OpCode, Site, Fishing_Grounds),
            by = "OpCode") %>%
  group_by(Site, Diet_Mouillot_2014) %>%
  summarise(n = n(), .groups="drop_last") %>%
  mutate(rel_abund = n / sum(n)) %>%
  pivot_wider(id_cols = Site,
              names_from = Diet_Mouillot_2014,
              values_from = rel_abund,
              values_fill = 0)

predictors <- left_join(fd_site, func_groups, by = "Site")

###############################################################################
# 2. Helper: convex hull vertices
###############################################################################

convhull.vert <- function(df) df[chull(df), , drop=FALSE]

###############################################################################
# 3. Build per-site convex hull plots
###############################################################################

titles  <- unique(fish_assemblage$Site)
palette <- viridis(length(titles), option="C")
t.space <- vector("list", length(titles))

hull.v  <- convhull.vert(point[,2:3])
hull.v2 <- convhull.vert(point[,4:5])

for (i in seq_along(titles)) {
  
  site_name <- titles[i]
  site_fg <- unique(fish_assemblage$Fishing_Grounds[fish_assemblage$Site == site_name])
  
  pointS <- fish_assemblage %>%
    filter(Site == site_name) %>%
    left_join(point, by="Species_full") %>%
    mutate(n = n / site.total$site.total[site.total$Site == site_name])
  
  if (nrow(pointS) == 0) {
    message("Skipping site with no species: ", site_name)
    next
  }
  
  s.hull.v  <- convhull.vert(pointS[, c("A1","A2")])
  s.hull.v2 <- convhull.vert(pointS[, c("A3","A4")])
  
  # Top 3 species abbreviated
  top_species <- pointS %>% arrange(desc(n)) %>% slice_head(n = 3)
  shortsp <- sapply(strsplit(top_species$Species_full, " "),
                    function(x) paste0(substr(x[1],1,1), ". ", x[2]))
  pointS$ShortSp <- ""
  pointS$ShortSp[match(top_species$Species_full, pointS$Species_full)] <- shortsp
  
  p1 <- ggplot() + theme_clean +
    geom_polygon(data=hull.v, aes(x=A1, y=A2), color='grey70', fill='white') +
    geom_polygon(data=s.hull.v, aes(x=A1, y=A2), fill=palette[i], alpha=0.4) +
    geom_point(data=pointS, aes(x=A1, y=A2, size=n),
               shape=21, fill=palette[i], color="grey30") +
    geom_text_repel(data=top_species, aes(x=A1, y=A2, label=shortsp),
                    size=2.3, fontface="italic", box.padding=0.5) +
    labs(x="PCo1", y="PCo2",
         title=paste0(site_name, "\n(", site_fg, ")")) +
    scale_size_continuous(range=c(1.5,8), guide="none")
  
  p2 <- ggplot() + theme_clean +
    geom_polygon(data=hull.v2, aes(x=A3, y=A4), color='grey70', fill='white') +
    geom_polygon(data=s.hull.v2, aes(x=A3, y=A4), fill=palette[i], alpha=0.4) +
    geom_point(data=pointS, aes(x=A3, y=A4, size=n),
               shape=21, fill=palette[i], color="grey30") +
    geom_text_repel(data=top_species, aes(x=A3, y=A4, label=shortsp),
                    size=2.3, fontface="italic", box.padding=0.5) +
    labs(x="PCo3", y="PCo4") +
    scale_size_continuous(range=c(1.5,8), guide="none")
  
  t.space[[i]] <- p1 / p2
}

###############################################################################
# 4. Global trait space (all sites combined)
###############################################################################

Sp.count <- fish_assemblage %>%
  group_by(Species_full) %>%
  summarise(n = sum(n), .groups="drop")

point$n <- Sp.count$n[match(point$Species_full, Sp.count$Species_full)]
point$n <- point$n / sum(point$n, na.rm=TRUE)

global1 <- ggplot() + theme_clean +
  geom_polygon(data=hull.v, aes(x=A1, y=A2), fill="transparent", color="grey") +
  geom_point(data=left_join(point, traits, by="Species_full"),
             aes(x=A1, y=A2, shape=ForageMode), size=2, color="black", alpha=0.6) +
  labs(x="PCo1", y="PCo2", title="All sites") +
  scale_shape_manual(values=1:10, name="Foraging mode")

global2 <- ggplot() + theme_clean +
  geom_polygon(data=hull.v2, aes(x=A3, y=A4), fill="transparent", color="grey") +
  geom_point(data=left_join(point, traits, by="Species_full"),
             aes(x=A3, y=A4, shape=ForageMode), size=2, color="black", alpha=0.6) +
  labs(x="PCo3", y="PCo4") +
  scale_shape_manual(values=1:10, guide="none")

###############################################################################
# 5. Barplots for functional diversity and diet guilds
###############################################################################

# Keep only relevant columns (drop OpCode etc.)
pred_long <- predictors %>%
  select(Site, Fishing_Grounds, TOP, TEve, TDiv,
         FC, HD, HM, IM, IS, OM, PK) %>%  # adjust to your actual trait/diet codes present
  pivot_longer(cols = -c(Site, Fishing_Grounds),
               names_to = "predictor",
               values_to = "value")

FD_bar <- pred_long %>%
  filter(predictor %in% c("TOP","TEve","TDiv")) %>%
  ggplot(aes(x=predictor, y=value, fill=Site)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  labs(y="Index measure", x=NULL) +
  scale_fill_viridis_d(option="C", guide="none") +
  theme_clean +
  scale_y_continuous(expand=expansion(mult=c(0,0.05)), limits=c(0,1))

A_bar <- pred_long %>%
  filter(!predictor %in% c("TOP","TEve","TDiv")) %>%
  ggplot(aes(x=predictor, y=value, fill=Site)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  labs(y="Relative abundance", x=NULL) +
  scale_fill_viridis_d(option="C", guide="none") +
  theme_clean +
  scale_y_continuous(expand=expansion(mult=c(0,0.05)))

bars <- (A_bar | FD_bar) * plot_layout(widths=c(2,3))

###############################################################################
# 6. Order sites by Fishing Grounds and combine figure
###############################################################################

site_order <- opcode_groups %>%
  distinct(Site, 'Fishing Grounds') %>%
  arrange('Fishing Grounds') %>%
  pull(Site)

t.space <- t.space[match(site_order, titles)]

sites_panel <- wrap_plots(t.space, ncol = 4)

Fig2 <- (bars /
           (global1 / global2) /
           sites_panel) *
  plot_layout(guides = "collect",
              heights = c(1.3, 2, 8)) &
  theme(plot.title = element_text(face = "bold",
                                  hjust = 0.5,
                                  size = 11))

Fig2

###############################################################################
# 7. Export (A3 portrait, thesis quality)
###############################################################################

ggsave("TraitSpaces_Thesis.svg",
       plot = Fig2,
       width = 42, height = 60, units = "cm", dpi = 300)
###############################################################################