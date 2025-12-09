library(mFD)

#prep fish DOV dataset
clean_data <- read.csv("Fish_Data_Fiji_Clean.csv")
clean_data$Transect <- as.factor(clean_data$Transect)

dataAll<- read.csv("DOV_All_EM_biomass.csv")
dataAll$Period <- as.factor(dataAll$Period)
head(dataAll)

common_levels <- intersect(levels(dataAll$Period), levels(clean_data$Transect))

dataAll$Period <- factor(dataAll$Period, levels = common_levels)
clean_data$Transect <- factor(clean_data$Transect, levels = common_levels)
dataAllF <- dataAll[dataAll$Period %in% clean_data$Transect, ]

dataAllF <- dataAllF %>%
  left_join(
    clean_data %>% select(Transect, Fishing_Ground),
    by = c("Period" = "Transect")
  )

#prep_matrix for fd
bio_matrix <- dataAllF %>%
  group_by(Period, Species_full) %>%
  summarise(biomass = sum(biomass), .groups = "drop") %>%
  pivot_wider(names_from = Species_full, values_from = biomass, values_fill = 0) %>%
  column_to_rownames("Period")

#prep_matrix fat fishing ground level for plot purpouses
FG_matrix <- dataAllF %>%
  group_by(Fishing_Ground, Species_full) %>%
  summarise(biomass = sum(biomass), .groups = "drop") %>%
  pivot_wider(names_from = Species_full, values_from = biomass, values_fill = 0) %>%
  column_to_rownames("Fishing_Ground")


#Trait databse (Created by Peter)
traits_mat<-read.csv("Fiji_traits.csv")

#matching datasets
traits_mat_unique <- traits_mat[!duplicated(traits_mat), ]
rownames(traits_mat_unique) <- traits_mat_unique$X
traits_mat_unique<-traits_mat_unique[,-1]
dim(traits_mat_unique)

common_species <- intersect(rownames(traits_mat_unique), colnames(bio_matrix))
traits_mat2 <- traits_mat_unique[common_species, ]
bio_matrix2 <-bio_matrix[, common_species]
FG_matrix2<--FG_matrix[, common_species]

# Convert categorical traits to factors
traits_mat2 <- traits_mat2 %>%
  mutate(
    Mobility = as.ordered(Mobility),
    Activity = as.ordered(Activity),
    Schooling = as.ordered(Schooling),
    Position = as.ordered(Position),
    Diet_Mouillot_2014 = as.factor(Diet_Mouillot_2014)
  )

#a data frame summarizing species traits (species in rows, traits in columns). 
fish_traits <- traits_mat2

#a matrix summarizing species gathering into assemblages (assemblages in rows, species in columns).
baskets_fish_weightsB<- as.matrix(bio_matrix2)
baskets_fish_weightsFG<- as.matrix(FG_matrix2)

#identical(colnames(baskets_fish_weightsFG),rownames(fish_traits))

#a data frame summarizing traits category (first column with traits name, second column with traits type, 
#third column with fuzzy name of fuzzy traits - if no fuzzy traits: NA).
fish_traits_cat <- read.csv("meta_trait.csv")

fish_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = fish_traits_cat,   
  sp_tr      = fish_traits, 
  stop_if_NA = TRUE)


sp_dist_fish <- mFD::funct.dist(
  sp_tr         = fish_traits,
  tr_cat        = fish_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

fspaces_quality_fish <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fish,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

round(fspaces_quality_fish$"quality_fspaces", 3)            # Quality metrics of spaces


mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fish,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d",
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")


sp_faxes_coord_fish <- fspaces_quality_fish$"details_fspaces"$"sp_pc_coord"

fish_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = fish_traits, 
  sp_faxes_coord = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)


# Print traits with significant effect:
fish_tr_faxes$"tr_faxes_stat"[which(fish_tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]


sp_faxes_coord_fish <- fspaces_quality_fish$"details_fspaces"$"sp_pc_coord"



#Calculate fd metrics for LMMs 
alpha_fd_indices_fish <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_fish_weightsB,
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_fish <- alpha_fd_indices_fish$"functional_diversity_indices"
fd_ind_values_fish$Transect <- rownames(fd_ind_values_fish)


clean_data <- clean_data %>%
  left_join(
    fd_ind_values_fish, by = c("Transect")
  )

#write.csv(clean_data, "Fish_Data_Fiji_Final_LF_Dez2024.csv")

#GO TO MODELS SCRIPT - RUN THE MODELS

#PLOTS

library(ggplot2)
library(ggrepel)

# Step 1: Extract coordinates and species names
coords <- as.data.frame(sp_faxes_coord_fish[, c("PC3", "PC4")])
coords$Species <- rownames(coords)

# Step 2: Find convex hull indices and order them correctly
hull_indices <- chull(coords$PC3, coords$PC4)
hull_indices <- c(hull_indices, hull_indices[1])  # close the polygon

# Step 3: Create a clean data frame for the hull
hull_coords <- coords[hull_indices, ]
hull_coords <- as.data.frame(hull_coords)  # ensure clean data frame
rownames(hull_coords) <- NULL

# Step 4: Plot
p <- ggplot(coords, aes(x = PC3, y = PC4)) +
  # All points
  geom_point(color = "darkgreen", fill = "white", size = 2, shape = 21) +
  # Convex hull polygon
  geom_polygon(data = hull_coords, aes(x = PC3, y = PC4),
               fill = "white", color = "black", alpha = 0.5) +
  # Hull vertices
  geom_point(data = hull_coords, aes(x = PC3, y = PC4),
             shape = 23, size = 3, fill = "blueviolet", color = "blueviolet") +
  # Vertex labels
  geom_text_repel(data = hull_coords, aes(label = Species),
                  size = 3, color = "black", fontface = "plain",
                  max.overlaps = Inf) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "grey95", color = NA))

p + theme_minimal()


View(coords)

# Step 5: Display plot




