# Run after 2) Functional Space
# SIZE FREQUENCY DISTRIBUTION
# Version date - 17/11/25
plot_size_frequency <- function(data, 
                                length_col = "Length_cm", 
                                group_col = NULL, 
                                binwidth = 2, 
                                relative = FALSE) {
  if (!length_col %in% names(data)) stop(paste("Column", length_col, "not found"))
  if (!is.null(group_col) && !group_col %in% names(data)) stop(paste("Grouping column", group_col, "not found"))
  
  # Summarise means/medians
  if (!is.null(group_col)) {
    stats_df <- data %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        mean_len = mean(.data[[length_col]], na.rm = TRUE),
        median_len = median(.data[[length_col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(!is.na(mean_len), !is.na(median_len))
  } else {
    stats_df <- tibble(
      mean_len = mean(data[[length_col]], na.rm = TRUE),
      median_len = median(data[[length_col]], na.rm = TRUE)
    )
  }
  
  p <- ggplot(data, aes(x = .data[[length_col]])) +
    geom_histogram(
      aes(y = if (relative) after_stat(count / sum(count)) else after_stat(count)),
      binwidth = binwidth,
      boundary = 0,
      fill = "skyblue",
      color = "grey30",
      alpha = 0.8
    ) +
    geom_vline(data = stats_df, aes(xintercept = mean_len), linetype = "dashed", color = "red", size = 0.9) +
    geom_vline(data = stats_df, aes(xintercept = median_len), linetype = "dotted", color = "blue", size = 0.9) +
    labs(
      x = "Total Length (cm)",
      y = if (relative) "Relative Frequency" else "Frequency",
      title = if (is.null(group_col)) "Size–Frequency Distribution" else paste("Size–Frequency by", group_col),
      subtitle = "Red dashed = Mean | Blue dotted = Median"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(group_col)) {
    p <- p + facet_wrap(as.formula(paste("~", group_col)), scales = if (relative) "fixed" else "free_y")
  }
  
  return(p)
}

# Overall size–frequency distribution - community-wide length structure.
p_overall <- plot_size_frequency(
  EmobsCombined_Filled,
  length_col = "Length_cm",
  binwidth = 2
) +
  labs(
    title = "Overall Size–Frequency Distribution of Reef Fishes",
    caption = "Red dashed = Mean length | Blue dotted = Median length"
  )

p_overall

# Filter to key target families - ecological differences in body size structure among key taxa.
target_families <- c(
  "Carangidae", "Caesionidae", "Labridae", "Lethrinidae", "Nemipteridae",
  "Siganidae", "Lutjanidae", "Mullidae", "Haemulidae", "Bodianinae",
  "Epinephelinae", "Scaridae", "Acanthuridae", "Leiognathidae"
)

fish_subset <- EmobsCombined_Filled %>%
  filter(Family %in% target_families) %>%
  mutate(Length_cm = as.numeric(Length_cm))

# Plot by family
p_family <- plot_size_frequency(
  fish_subset,
  length_col = "Length_cm",
  group_col = "Family",
  binwidth = 2
) +
  labs(
    title = "Size–Frequency Distributions of Target Fish Families",
    subtitle = "Red dashed = Mean | Blue dotted = Median",
    caption = "Facets show family-specific size structures; note differing y-axes."
  )

p_family

ggsave("SizeFrequency_Overall.png", p_overall, width = 8, height = 5, dpi = 400)
ggsave("SizeFrequency_ByFamily.png", p_family, width = 10, height = 8, dpi = 400)