dataAll<- read.csv("DOV_All_EM_biomass.csv")
dataAll$Period <- as.factor(dataAll$Period)

clean_data<- read.csv("Fish_Data_Fiji_Final_LF_Dez2024.csv")
clean_data$Transect <- as.factor(clean_data$Transect)

common_levels <- intersect(levels(dataAll$Period), levels(clean_data$Transect))

dataAll$Period <- factor(dataAll$Period, levels = common_levels)
clean_data$Transect <- factor(clean_data$Transect, levels = common_levels)
dataAllF <- dataAll[dataAll$Period %in% clean_data$Transect, ]

dataAllF <- dataAllF %>%
  left_join(
    clean_data %>% select(Transect, Fishing_Ground),
    by = c("Period" = "Transect")
  )

dataAllF$Fishing_Ground<-as.factor(dataAllF$Fishing_Ground)
levels(dataAllF$Fishing_Ground)

## Biomass

#Group biomass by species per transect

ComStData<- aggregate(biomass ~ Species_full + Period + Fishing_Ground, data=dataAllF, sum)
head(ComStData)

comm_matrix <- ComStData %>%
  pivot_wider(
    names_from = Species_full,          # new columns = species
    values_from = biomass,       # cell values = abundance
    values_fill = list(biomass = 0)  # replace missing abundance with 0
  )

meta <- ComStData %>%
  distinct(Period, Fishing_Ground) %>%
  arrange(match(Period, comm_matrix$Period))

head(comm_matrix)
comm_matrix <- comm_matrix[order(comm_matrix$Period), ]
meta <- meta[order(meta$Period), ]
rownames(comm_matrix) <- comm_matrix$Period
comm_matrix <- comm_matrix[ , -c(1,2)]
dist_bc <- vegdist(comm_matrix, method = "bray")


permanova_result <- adonis2(dist_bc ~ Fishing_Ground, data = meta, permutations = 999)
#save as csv 

#write.csv(permanova_result, "permanova_FishCommunity_result.csv")


#pairwise_perm <- pairwise.adonis2(dist_bc, meta, permutations = 999, method = "bray", p_adjust = "BH")

# --- Prepare ---
# comm_matrix: rows = samples, columns = species, values = abundances
# meta: dataframe with Sample IDs and Fishing_Ground

groups <- unique(meta$Fishing_Ground)
results <- data.frame()

# --- Loop through all pairs of Fishing_Ground ---
for (i in 1:(length(groups)-1)) {
  for (j in (i+1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    
    # Subset samples for these two groups
    idx <- meta$Fishing_Ground %in% c(g1, g2)
    sub_comm <- comm_matrix[idx, ]
    sub_meta <- meta[idx, , drop = FALSE]
    
    # Calculate Bray–Curtis dissimilarity
    sub_dist <- vegdist(sub_comm, method = "bray")
    
    # Run PERMANOVA
    test <- adonis2(sub_dist ~ Fishing_Ground, data = sub_meta, permutations = 999)
    
    #pairwise_perm <- pairwise.adonis2(sub_dist, analysis_data, permutations = 999, method = "bray", p_adjust = "BH")
    
    # Store results
    results <- rbind(
      results,
      data.frame(
        Group1 = g1,
        Group2 = g2,
        F = test$F[1],
        R2 = test$R2[1],
        p = test$`Pr(>F)`[1]
      )
    )
  }
}

# --- Adjust p-values for multiple comparisons (optional) ---
results$p.adj <- p.adjust(results$p, method = "BH")  # Benjamini-Hochberg

# --- View results ---
print(results)

#write.csv(results, "permanova_FishCommunity_pairwise_result.csv")



#NMDS for visualisation
set.seed(123)
nmds <- metaMDS(comm_matrix, distance = "bray", k = 2, trymax = 100)

nmds_points <- as.data.frame(nmds$points)
nmds_points$Fishing_Ground <- meta$Fishing_Ground

nmdsRFB<- ggplot(nmds_points, aes(MDS1, MDS2, color = Fishing_Ground)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 0.8) +
  scale_color_viridis_d(option = "F") +   # viridis discrete palette, option F
  theme_minimal(base_size = 12) +
  labs(
    title =  "Reef Fish Assemblages",
    x = "NMDS1",
    y = "NMDS2",
    color = "Fishing Ground"
  )




#Funtional Diversity

comm_matrix <- clean_data[,c("Transect","fric","fdiv","fspe","fori","feve","fdis")]

meta <- clean_data %>%
  distinct(Transect, Fishing_Ground) %>%
  arrange(match(Transect, comm_matrix$Transect))

head(comm_matrix)
comm_matrix <- comm_matrix[order(comm_matrix$Transect), ]
meta <- meta[order(meta$Transect), ]
rownames(comm_matrix) <- comm_matrix$Transect
comm_matrix <- comm_matrix[ , -c(1,2)]
dist_bc <- vegdist(comm_matrix, method = "euclidian")


# --- Prepare ---
# comm_matrix: rows = samples, columns = species, values = abundances
# meta: dataframe with Sample IDs and Fishing_Ground

groups <- unique(meta$Fishing_Ground)
results <- data.frame()

# --- Loop through all pairs of Fishing_Ground ---
for (i in 1:(length(groups)-1)) {
  for (j in (i+1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    
    # Subset samples for these two groups
    idx <- meta$Fishing_Ground %in% c(g1, g2)
    sub_comm <- comm_matrix[idx, ]
    sub_meta <- meta[idx, , drop = FALSE]
    
    # Calculate Bray–Curtis dissimilarity
    sub_dist <- vegdist(sub_comm, method = "euclidian")
    
    # Run PERMANOVA
    test <- adonis2(sub_dist ~ Fishing_Ground, data = sub_meta, permutations = 999)
    
    # Store results
    results <- rbind(
      results,
      data.frame(
        Group1 = g1,
        Group2 = g2,
        F = test$F[1],
        R2 = test$R2[1],
        p = test$`Pr(>F)`[1]
      )
    )
  }
}

# --- Adjust p-values for multiple comparisons (optional) ---
results$p.adj <- p.adjust(results$p, method = "BH")  # Benjamini-Hochberg

# --- View results ---
print(results)



set.seed(123)
nmds <- metaMDS(comm_matrix, distance = "euclidian", k = 2, trymax = 100)

nmds_points <- as.data.frame(nmds$points)
nmds_points$Fishing_Ground <- meta$Fishing_Ground

nmdsFD<-ggplot(nmds_points, aes(MDS1, MDS2, color = Fishing_Ground)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 0.8) +
  scale_color_viridis_d(option = "F") +   # viridis discrete palette, option F
  theme_minimal(base_size = 12) +
  labs(
    title = "Functional Diversity",
    x = "NMDS1",
    y = "NMDS2",
    color = "Fishing Ground"
  )



get_legend <- function(p) {
  tmp <- ggplotGrob(p)
  leg <- gtable::gtable_filter(tmp, "guide-box")
  return(leg)
}

shared_legend <- get_legend(nmdsFD)


nmdsFD_noleg <- nmdsFD + theme(legend.position = "none")
nmdsFC_noleg <- nmdsRFB + theme(legend.position = "none")

grid.arrange(
  arrangeGrob(nmdsFC_noleg,nmdsFD_noleg, nrow = 2),
  shared_legend,
  ncol = 2,
  widths = c(6, 3)   # adjust legend width if needed
)
