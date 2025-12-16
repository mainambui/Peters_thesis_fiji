
# ============================================================
# Functional space panels (a–e), Mouillot-style
# Global hull + ground hulls; points sized by biomass (0 = absence)
# Colours: viridis option "F" (rocket) for hulls AND points
# Fish shapes in ALL panels; NO text labels
# Figure-level legend; all panels identical size
# ============================================================

# 1) Packages -------------------------------------------------
pkgs <- c("tidyverse", "ggplot2", "viridis", "cowplot", "fishualize", "scales")
invisible(lapply(pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)))
lapply(pkgs, library, character.only = TRUE)

# 2) Load data -----------------------------------------------
pca_file <- "species_pca.csv"    # rows = species; cols = PC1, PC2
fg_file  <- "fg_matrix.csv"      # rows = grounds; cols = species; biomass; 0 = absence
tran_file <- "bio_matrix.csv"

spp_space <- read.csv(pca_file, row.names = 1, check.names = FALSE) %>%
  rownames_to_column("species") %>%
  mutate(species = trimws(gsub('^\"|\"$', '', species))) %>%
  select(species, PC1, PC2)

fg_mat <- read.csv(fg_file, row.names = 1, check.names = FALSE)
grounds <- rownames(fg_mat)

biom_long <- fg_mat %>%
  as.data.frame() %>%
  rownames_to_column("ground") %>%
  pivot_longer(-ground, names_to = "species", values_to = "biomass") %>%
  mutate(species = trimws(gsub('^\"|\"$', '', species))) %>%
  left_join(spp_space, by = "species")

stopifnot(all(!is.na(biom_long$PC1)), all(!is.na(biom_long$PC2)))

bio_mat <- read.csv(tran_file, row.names = 1, check.names = FALSE)
bio_mat<-bio_mat/1000 
Transect <- rownames(bio_mat)

biom_long2 <- bio_mat %>%
  as.data.frame() %>%
  rownames_to_column("Transect") %>%
  pivot_longer(-Transect, names_to = "species", values_to = "biomass") %>%
  mutate(species = trimws(gsub('^\"|\"$', '', species))) %>%
  left_join(spp_space, by = "species")

stopifnot(all(!is.na(biom_long2$PC1)), all(!is.na(biom_long2$PC2)))

#create columns "ground" in biom_long2 for aestethics in ggplot
clean_data <- read.csv("Fish_Data_Fiji_Clean.csv")
clean_data$Transect <- as.factor(clean_data$Transect)
clean_data<- clean_data[,c("Transect","Fishing_Ground")]

biom_long2 <-biom_long2 %>%
  as.data.frame() %>%
  left_join(
    clean_data,
    by = c("Transect")
  )

colnames(biom_long2)[6] <- "ground"

# 3) Convex hulls --------------------------------------------
get_hull <- function(df) {
  df <- dplyr::distinct(df, PC1, PC2, .keep_all = TRUE)
  if (nrow(df) < 3) return(df[0, ])
  df[chull(df$PC1, df$PC2), , drop = FALSE]
}

global_hull <- get_hull(spp_space)

ground_hulls <- biom_long %>%
  filter(biomass > 0) %>%                  # 0 = absence
  group_by(ground) %>%
  group_modify(~ get_hull(.x)) %>%
  ungroup()



# 4) Scales & palette ----------------------------------------
ground_levels <- grounds
biom_long    <- biom_long2    %>% mutate(ground = factor(ground, levels = ground_levels))
ground_hulls <- ground_hulls %>% mutate(ground = factor(ground, levels = ground_levels))


global_max_biomass <- max(biom_long$biomass, na.rm = TRUE)
size_breaks <- scales::pretty_breaks(n = 4)(c(0, global_max_biomass)) %>% unique() %>% round(1)

xlim_all <- range(spp_space$PC1) + c(-0.05, 0.05)
ylim_all <- range(spp_space$PC2) + c(-0.05, 0.05)

# Figure-wide size scale (identical across panels; legend extracted later)
size_scale_common <- scale_size_area(
  max_size = 7,
  limits   = c(0, global_max_biomass),
  breaks   = size_breaks,
  labels   = scales::label_number(),
  guide    = "none"   # per-panel guides hidden
)

# 5) Base layer ----------------------------------------------
base_layer <- ggplot() +
  geom_point(data = spp_space, aes(PC1, PC2), color = "grey40", size = 1.6, alpha = 0.1) +
  geom_polygon(data = global_hull, aes(PC1, PC2), fill = NA, color = "grey30", linewidth = 0.5) +
  coord_equal(xlim = xlim_all, ylim = ylim_all) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  labs(x = "PC1", y = "PC2")

# 6) Ground panels (legend hidden; colours matched) ----------
plots <- vector("list", length(ground_levels))
for (i in seq_along(ground_levels)) {
  g <- ground_levels[i]
  points_g <- biom_long %>% filter(ground == g, biomass > 0)
  hull_g   <- ground_hulls %>% filter(ground == g)
  
  p <- base_layer +
    geom_polygon(data = hull_g, aes(PC1, PC2, fill = ground), color = NA, alpha = 0.35) +
    geom_point(data = points_g, aes(PC1, PC2, size = biomass,fill = ground), color="black", shape=21, alpha = 0.92) +
    size_scale_common +
    scale_fill_viridis_d(option = "F", limits = ground_levels, drop = FALSE, guide = "none") +  # rocket
    scale_colour_viridis_d(option = "F", limits = ground_levels, drop = FALSE, guide = "none") +
    labs(title = as.character(g))
  
  plots[[i]] <- p
}

# 7) Silhouette policy (clean, deterministic, NO text) -------
# Species overrides (exact surrogates)
species_override_map <- c(
  "Plectropomus laevis" = "Cephalopholis_argus",   # Serranidae
  "Hemigymnus fasciatus" = "Epibulus_insidiator"   # Labridae
)

# Genus overrides (preferred silhouettes)
genus_override_map <- c(
  "Naso"        = "Naso_unicornis",               # true acanthurid
  "Halichoeres" = "Coris_gaimard",                # force Coris for Halichoeres sp
  "Chromis"     = "Acanthochromis_polyacanthus",  # Pomacentridae surrogate (no Chromis in catalog)
  "Pomacentrus" = "Acanthochromis_polyacanthus",
  "Gymnocaesio" = NA_character_                   # require Caesionidae; if none → no shape
)

# Genus → expected family (prevents wrong-family fallbacks)
genus_to_family <- c(
  "Naso"          = "Acanthuridae",
  "Halichoeres"   = "Labridae",
  "Hemigymnus"    = "Labridae",
  "Coris"         = "Labridae",
  "Chromis"       = "Pomacentridae",
  "Pomacentrus"   = "Pomacentridae",
  "Chrysiptera"   = "Pomacentridae",
  "Neopomacentrus"= "Pomacentridae",
  "Plectropomus"  = "Serranidae",
  "Monotaxis"     = "Lethrinidae",
  "Carangoides"   = "Carangidae",
  "Gymnocaesio"   = "Caesionidae",   # strictly Caesionidae
  "Taeniamia"     = "Apogonidae",
  "Cephalopholis" = "Serranidae",
  "Chaetodon"     = "Chaetodontidae",
  "Parupeneus"    = "Mullidae",
  "Kyphosus"      = "Kyphosidae",
  "Sargocentron"  = "Holocentridae",
  "Myripristis"   = "Holocentridae",
  "Arothron"      = "Tetraodontidae",
  "Zanclus"       = "Zanclidae",
  "Siganus"       = "Siganidae"
)

# Cache silhouettes catalog (family + option + genus)
sil_tbl <- fishualize::fishapes() %>%
  dplyr::mutate(genus_opt = sub("_.*$", "", option))   # "Coris_gaimard" -> "Coris")

# NA-safe equality helper
safe_is_expected <- function(fam, expect_family) {
  isTRUE(is.na(expect_family)) || isTRUE(!is.na(fam) && fam == expect_family)
}

# Helpers -----------------------------------------------------
resolve_option_to_family <- function(option_string) {
  row <- sil_tbl %>% dplyr::filter(option == option_string)
  if (nrow(row) > 0) row$family[1] else NA_character_
}
pick_first_in_genus <- function(genus) {
  if (is.na(genus)) return(list(family = NA_character_, option = NA_character_))
  rows <- sil_tbl %>% dplyr::filter(genus_opt == genus)
  if (nrow(rows) > 0) list(family = rows$family[1], option = rows$option[1]) else list(family = NA_character_, option = NA_character_)
}
pick_first_in_family <- function(family_name) {
  if (is.na(family_name)) return(list(family = NA_character_, option = NA_character_))
  rows <- sil_tbl %>% dplyr::filter(family == family_name)
  if (nrow(rows) > 0) list(family = rows$family[1], option = rows$option[1]) else list(family = NA_character_, option = NA_character_)
}

# Master resolver (returns silhouette or none; NO text) -------
get_fish_silhouette <- function(species_name) {
  sp     <- trimws(species_name)
  sp_key <- gsub(" +", "_", sp)
  genus  <- sub(" .*$", "", sp)
  expect_family <- if (genus %in% names(genus_to_family)) genus_to_family[[genus]] else NA_character_
  
  # 0) Species override
  if (sp %in% names(species_override_map)) {
    opt <- species_override_map[[sp]]
    fam <- resolve_option_to_family(opt)
    if (safe_is_expected(fam, expect_family)) {
      return(list(family = fam, option = opt, match_type = "species_override"))
    }
  }
  
  # 1) Exact species silhouette
  exact <- sil_tbl %>% dplyr::filter(option == sp_key)
  if (nrow(exact) > 0) {
    fam <- exact$family[1]
    if (safe_is_expected(fam, expect_family)) {
      return(list(family = fam, option = exact$option[1], match_type = "species"))
    }
  }
  
  # 2) Genus override (preferred)
  if (genus %in% names(genus_override_map)) {
    opt <- genus_override_map[[genus]]
    if (!is.na(opt)) {
      fam <- resolve_option_to_family(opt)
      if (safe_is_expected(fam, expect_family)) {
        return(list(family = fam, option = opt, match_type = "genus_override"))
      }
    } else {
      # NA means "require expected family only" (e.g., Gymnocaesio → Caesionidae only)
      fam_fb <- pick_first_in_family(expect_family)
      if (!is.na(fam_fb$option)) {
        return(list(family = fam_fb$family, option = fam_fb$option, match_type = "family"))
      } else {
        return(list(family = NA_character_, option = NA_character_, match_type = "none"))
      }
    }
  }
  
  # 3) Genus fallback (first available silhouette in that genus)
  gen_fb <- pick_first_in_genus(genus)
  if (!is.na(gen_fb$option) && safe_is_expected(gen_fb$family, expect_family)) {
    return(list(family = gen_fb$family, option = gen_fb$option, match_type = "genus"))
  }
  
  # 4) Family fallback (only within expected family)
  fam_fb <- pick_first_in_family(expect_family)
  if (!is.na(fam_fb$option)) {
    return(list(family = fam_fb$family, option = fam_fb$option, match_type = "family"))
  }
  
  # 5) None → no shape
  list(family = NA_character_, option = NA_character_, match_type = "none")
}

# Reporter: table of matches (saved for your records)
check_global_hull_silhouettes <- function(global_hull) {
  purrr::pmap_dfr(global_hull, function(species, PC1, PC2, ...) {
    m <- get_fish_silhouette(species)
    tibble(species = species, PC1 = PC1, PC2 = PC2,
           family = m$family, option = m$option, match_type = m$match_type)
  })
}

# Drawer: add ONLY silhouettes for each vertex (NO text), to ONE plot --------
add_global_hull_shapes <- function(plot_obj, hull_df,
                                   xlim_all = NULL, ylim_all = NULL,
                                   fill = "#2B2B2B", alpha = 0.6) {
  
  if (is.null(xlim_all)) xlim_all <- range(hull_df$PC1)
  if (is.null(ylim_all)) ylim_all <- range(hull_df$PC2)
  box_w <- diff(xlim_all) * 0.06
  box_h <- diff(ylim_all) * 0.09
  
  matches <- check_global_hull_silhouettes(hull_df)
  # Log once (only for the first call)
  if (!file.exists("global_hull_silhouettes_used.csv")) {
    readr::write_csv(matches, "global_hull_silhouettes_used.csv")
  }
  
  p <- plot_obj
  for (i in seq_len(nrow(matches))) {
    row <- matches[i, ]
    # Draw silhouette if available
    if (row$match_type %in% c("species_override","species","genus_override","genus","family") &&
        !is.na(row$family) && !is.na(row$option)) {
      p <- p + fishualize::add_fishape(
        family = row$family, option = row$option,
        xmin = row$PC1 - box_w/1, xmax = row$PC1 + box_w/1,
        ymin = row$PC2 - box_h/1, ymax = row$PC2 + box_h/1,
        fill = fill, alpha = alpha
      )
    }
  }
  p
}

# Apply silhouettes to ALL panels (a–e) ----------------------
for (i in seq_along(plots)) {
  plots[[i]] <- add_global_hull_shapes(plots[[i]], global_hull,
                                       xlim_all = xlim_all, ylim_all = ylim_all,
                                       fill = "#2B2B2B", alpha = 0.5)
}

# 8) Figure-level biomass legend (extract once, place below) --
legend_dummy <- ggplot(
  biom_long %>% filter(biomass > 0) %>% slice_head(n = 50),
  aes(x = PC1, y = PC2, size = biomass)
) +
  geom_point(alpha = 0.3, color="black",fill = "gray",shape=21) +                     # invisible points; legend only
  scale_size_area(
    max_size = 7,
    limits   = c(2, global_max_biomass),
    breaks   = size_breaks,
    labels   = scales::label_number(),
    guide    = guide_legend(title = "Biomass (Kg/250m2)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(5, 5, 5, 5),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9)
  )

legend <- cowplot::get_legend(legend_dummy)


# 9) Assemble panels (a–e) identically sized + legend ----------
panel_grid <- cowplot::plot_grid(plotlist = plots, labels = letters[1:5], ncol = 3)

final_fig <- cowplot::plot_grid(
  panel_grid,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.12)      # slightly more space to avoid legend clipping
)

ggsave("functional_space_panels.png", final_fig, width = 11.5, height = 8.7, dpi = 300)
print(final_fig)

# ============================================================
# End of script
