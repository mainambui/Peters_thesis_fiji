#read top models


modBio<- all_responses_lmer_logBiomass
modAbun <- all_responses_lmer_logabun
modFric<-all_responses_lmer_fric
modFeve<-all_responses_lm_feve

library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)

library(MuMIn)
library(performance)

r2(all_responses_lmer_logBiomass)
r2(all_responses_lmer_logabun)
r2(all_responses_lmer_fric)


#-----------------------------
# Function to prepare plotting data from lm
#-----------------------------
prepare_plot_data <- function(lm_model, model_name) {
  coef_df <- as.data.frame(coef(summary(lm_model)))
  coef_df <- coef_df %>%
    rownames_to_column(var = "term") %>%
    rename(estimate = "Estimate",
           std_error = "Std. Error",
           t_value = "t value",
           p_value = "Pr(>|t|)") %>%
    # Exclude intercept
    filter(term != "(Intercept)") %>%
    mutate(
      conf.low  = estimate - 1.96 * std_error,
      conf.high = estimate + 1.96 * std_error,
      sig_color = ifelse(p_value < 0.05, "darkgreen", "lightgrey"),
      model = model_name
    )
  
  return(coef_df)
}

#-----------------------------
# Prepare data for all models
#-----------------------------
dfBio   <- prepare_plot_data(modBio, "Fish Biomass")
dfAbun  <- prepare_plot_data(modAbun, "Fish Abundance")
dfFric  <- prepare_plot_data(modFric, "Functional Richness")
dfFeve  <- prepare_plot_data(modFeve, "Funcional Evenness")

# Combine all into one data frame
df_all <- bind_rows(dfBio, dfAbun, dfFric, dfFeve)

#-----------------------------
# Plot side by side
#-----------------------------
plot_single_model <- function(df_model) {
  ggplot(df_model, aes(x = estimate, y = term, color = sig_color)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    scale_color_identity() +
    theme_bw() +
    labs(
      x = "Estimate",
      y = "Predictor",
      title = paste0("", unique(df_model$model))
    )
}



pBio   <- plot_single_model(dfBio)
pAbun  <- plot_single_model(dfAbun)
pFric  <- plot_single_model(dfFric)
pFeve  <- plot_single_model(dfFeve)

#-----------------------------
# Display plots
#-----------------------------
library(gridExtra)
grid.arrange(pBio, pAbun, pFric, pFeve,
             ncol = 2, nrow = 2)


#####

library(tidyverse)
library(viridis)

df_long <- clean_data %>%
  select(Fishing_Ground,
         biomass_kg, abundance_ind_250m2) %>%
  pivot_longer(
    cols = -Fishing_Ground,
    names_to = "Variable",
    values_to = "Value"
  )


ggplot(df_long, aes(x = Fishing_Ground, y = Value, fill = Fishing_Ground)) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~ Variable, scales = "free_y") +
  scale_fill_viridis_d(option = "F") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),      # ← remove x-axis labels
    axis.ticks.x = element_blank(),     # ← optional: remove ticks too
    legend.position = "right"
  ) +
  labs(
    x = "",
    y = "",
    fill = "Fishing Ground"
  ) + theme_classic()

