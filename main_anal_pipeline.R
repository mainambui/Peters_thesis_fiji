#  ANALYSIS PIPELINE 

# Date: 2025
# Description: Complete analysis pipeline for reef ecology data including:
#              - Data validation and cleaning
#              - Linear and mixed effects models
#              - Multivariate analysis (PCA, PERMANOVA)
#              - Non-parametric tests
#              - Comprehensive visualization and reporting

# OUTPUTS:
# - figures/: Visualizations and plots
# - tables/: Statistical results and coefficients
# - models/: Model objects and summaries
# - ordination/: Multivariate analysis results
# - diagnostics/: Model diagnostic plots
# - reports/: Analysis summaries

# CLEAR WORKSPACE AND LOAD LIBRARIES
# ==================================
rm(list = ls())
cat("Initializing Reef Ecology Analysis Pipeline...\n")

# Function to install and load required packages
install_and_load_packages <- function() {
  required_packages <- c(
    "tidyverse", "lme4", "lmerTest", "vegan", "ggplot2",
    "car", "MuMIn", "DHARMa", "ggpubr", "GGally",
    "broom", "broom.mixed", "glue", "ggcorrplot"
  )
  
  # Install missing packages
  inst <- rownames(installed.packages())
  missing_packages <- setdiff(required_packages, inst)
  if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, dependencies = TRUE)
  }
  
  # Load all packages
  cat("Loading required packages...\n")
  suppressPackageStartupMessages({
    library(tidyverse)
    library(lme4)
    library(lmerTest)
    library(vegan)
    library(ggplot2)
    library(car)
    library(MuMIn)
    library(DHARMa)
    library(ggpubr)
    library(GGally)
    library(broom)
    library(broom.mixed)
    library(glue)
    library(ggcorrplot)
  })
  
  cat("✅ All packages loaded successfully!\n")
}

#Install and load packages
install_and_load_packages()

#SET WORKING DIRECTORY AND LOAD DATA

cat("Setting up working directory and loading data...\n")

# Check if working directory is set, if not use current (user-editable default)
if (!exists("custom_wd")) {
  custom_wd <- getwd()
  cat("Using working directory:", custom_wd, "\n")
}
if (!dir.exists(custom_wd)) {
  stop("Working directory does not exist: ", custom_wd)
}
setwd(custom_wd)

# Load data with error handling
load_reef_data <- function(file_path = "DOV_transect_connect_final.csv") {
  if (file.exists(file_path)) {
    cat("Loading data from:", file_path, "\n")
    data <- read.csv(file_path)
    cat("ata loaded successfully:", nrow(data), "rows,", ncol(data), "columns\n")
    return(data)
  } else {
    stop("Data file not found: ", file_path,
         "\nPlease ensure the file exists in the working directory.")
  }
}

data <- load_reef_data()

# ENHANCED DATA VALIDATION FUNCTION
# =================================
validate_analysis_data <- function(data, response_vars, predictor_vars) {
  cat("=== DATA VALIDATION CHECKS ===\n")
  issues <- list(warnings = character(), errors = character())
  
  # Check data existence
  if (nrow(data) == 0) {
    issues$errors <- c(issues$errors, "Dataset is empty")
    return(issues)
  }
  
  # Check variable existence
  missing_vars <- c(
    setdiff(response_vars, names(data)),
    setdiff(predictor_vars, names(data))
  )
  
  if (length(missing_vars) > 0) {
    issues$errors <- c(issues$errors, paste("Missing variables:", paste(missing_vars, collapse = ", ")))
  }
  
  # Check for sufficient sample size
  if (nrow(data) < 10) {
    issues$warnings <- c(issues$warnings, "Warning: Small sample size (n < 10) may affect model reliability")
  }
  
  # Check variance in response variables + excessive zeros
  for (var in response_vars) {
    if (var %in% names(data)) {
      var_vector <- data[[var]]
      if (all(is.na(var_vector))) {
        issues$warnings <- c(issues$warnings, paste("All values are NA in response variable:", var))
      } else {
        if (var(var_vector, na.rm = TRUE) == 0) {
          issues$warnings <- c(issues$warnings, paste("Zero variance in response variable:", var))
        }
        denom <- sum(!is.na(var_vector))
        if (denom > 0) {
          zero_prop <- sum(var_vector == 0, na.rm = TRUE) / denom
          if (zero_prop > 0.8) {
            issues$warnings <- c(issues$warnings, paste0("High proportion of zeros (", round(zero_prop * 100, 1), "%) in: ", var))
          }
        }
      }
    }
  }
  
  # Check factor levels (guard for non-existing vars)
  existing_preds <- predictor_vars[predictor_vars %in% names(data)]
  factor_flags <- sapply(existing_preds, function(x) is.factor(data[[x]]))
  factor_vars <- existing_preds[factor_flags]
  for (fac in factor_vars) {
    if (length(unique(na.omit(data[[fac]]))) < 2) {
      issues$errors <- c(issues$errors, paste("Factor", fac, "has only one level"))
    }
  }
  
  #Check for complete cases (only over available columns)
  cc_cols <- c(intersect(response_vars, names(data)), intersect(predictor_vars, names(data)))
  if (length(cc_cols) > 0) {
    complete_cases <- complete.cases(data[, cc_cols, drop = FALSE])
    if (sum(complete_cases) < nrow(data) * 0.5) {
      issues$warnings <- c(
        issues$warnings,
        paste("Only", sum(complete_cases), "complete cases out of", nrow(data), "observations")
      )
    }
  }
  
  #Print validation results
  if (length(issues$errors) > 0) {
    cat("RRORS:\n"); cat(paste("  -", issues$errors), sep = "\n")
  }
  if (length(issues$warnings) > 0) {
    cat(" WARNINGS:\n"); cat(paste("  -", issues$warnings), sep = "\n")
  }
  if (length(issues$errors) == 0 & length(issues$warnings) == 0) {
    cat("All data validation checks passed\n")
  }
  
  return(issues)
}

#DATA CLEANING AND PREPARATION

cat("=== DATA CLEANING AND PREPARATION ===\n")
clean_data <- data %>%
  #Convert character columns to factors
  mutate(across(where(is.character), as.factor)) %>%
  #Remove rows with missing values in key columns (guard columns that may not exist)
  filter(!if_any(any_of(c("biomass_kg_ha", "abundance_ind_250m2", "FRic", "FEve", "FDiv", "FDis", "RaoQ")), is.na)) %>%
  #Create unique site identifier
  mutate(Site_ID = paste(Location, Dive_Site, sep = "_")) %>%
  #Keep only locations with sufficient replicates
  group_by(Location) %>%
  filter(n() >= 3) %>%
  ungroup()

cat("Data cleaning completed:\n")
cat("  - Original data:", nrow(data), "rows\n")
cat("  - Cleaned data:", nrow(clean_data), "rows\n")
cat("  - Locations:", paste(unique(clean_data$Location), collapse = ", "), "\n")

#CREATE OUTPUT DIRECTORIES

create_output_directories <- function() {
  dirs <- c("figures", "tables", "models", "ordination", "diagnostics", "reports")
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, showWarnings = FALSE, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
  cat("ll output directories created successfully!\n")
}
create_output_directories()

#ENHANCED PAIRWISE PERMANOVA FUNCTION

pairwise.adonis2 <- function(x, data, permutations = 999, method = "euclidean", p_adjust = "BH", ...) {
  locations <- unique(data$Location)
  results <- list()
  pairwise_results <- data.frame()
  
  for (i in 1:(length(locations) - 1)) {
    for (j in (i + 1):length(locations)) {
      loc1 <- locations[i]
      loc2 <- locations[j]
      
      #Subset data for this pair
      pair_data <- data[data$Location %in% c(loc1, loc2), , drop = FALSE]
      pair_matrix <- x[data$Location %in% c(loc1, loc2), , drop = FALSE]
      
      if (nrow(pair_matrix) > 2 && ncol(pair_matrix) > 0) {
        dist_matrix <- vegan::vegdist(pair_matrix, method = method)
        permanova <- vegan::adonis2(dist_matrix ~ Location, data = pair_data, permutations = permutations)
        
        result_row <- data.frame(
          Comparison = paste(loc1, "vs", loc2),
          F_value = permanova$F[1],
          R2 = permanova$R2[1],
          p_value = permanova$`Pr(>F)`[1],
          n = nrow(pair_data)
        )
        
        pairwise_results <- rbind(pairwise_results, result_row)
        results[[paste(loc1, "vs", loc2)]] <- list(
          F_value = permanova$F[1],
          R2 = permanova$R2[1],
          p_value = permanova$`Pr(>F)`[1],
          n = nrow(pair_data)
        )
      }
    }
  }
  
  #Apply p-value adjustment
  if (nrow(pairwise_results) > 0) {
    pairwise_results$p_adj <- p.adjust(pairwise_results$p_value, method = p_adjust)
    
    #Update results list with adjusted p-values
    for (i in 1:nrow(pairwise_results)) {
      comp_name <- pairwise_results$Comparison[i]
      results[[comp_name]]$p_adj <- pairwise_results$p_adj[i]
    }
  }
  
  return(list(
    pairwise_results = pairwise_results,
    detailed_results = results
  ))
}

#MODEL DIAGNOSTICS FUNCTION

create_model_diagnostics <- function(model, response_var, model_type = "lm") {
  plots <- list()
  diagnostics <- list()
  
  if (inherits(model, "lm") | inherits(model, "lmerMod")) {
    
    #Extract residuals and fitted values
    if (inherits(model, "lmerMod")) {
      fitted_vals <- fitted(model)
      residuals_vals <- residuals(model, type = "pearson")
    } else {
      fitted_vals <- fitted(model)
      residuals_vals <- residuals(model)
    }
    
    #1.Residuals vs Fitted plot
    p1 <- ggplot(data = data.frame(fitted = fitted_vals,
                                   residuals = residuals_vals),
                 aes(x = fitted, y = residuals)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "loess", se = TRUE, color = "blue", alpha = 0.2) +
      labs(title = paste("Residuals vs Fitted -", response_var),
           subtitle = paste("Model type:", model_type),
           x = "Fitted values", y = "Residuals") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    #2.Q-Q plot for normality
    p2 <- ggplot(data = data.frame(residuals = residuals_vals),
                 aes(sample = residuals)) +
      stat_qq(alpha = 0.6, size = 2) +
      stat_qq_line(color = "red") +
      labs(title = paste("Q-Q Plot -", response_var),
           subtitle = "Checking normality assumption",
           x = "Theoretical Quantiles", y = "Sample Quantiles") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    #3.Scale-Location plot
    p3 <- ggplot(data = data.frame(fitted = fitted_vals,
                                   sqrt_abs_resid = sqrt(abs(residuals_vals))),
                 aes(x = fitted, y = sqrt_abs_resid)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "loess", se = TRUE, color = "blue", alpha = 0.2) +
      labs(title = paste("Scale-Location Plot -", response_var),
           subtitle = "Checking homoscedasticity",
           x = "Fitted values", y = "√|Standardized residuals|") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    #4.Residuals vs Leverage
    if (inherits(model, "lm")) {
      leverage <- hatvalues(model)
      p4 <- ggplot(data = data.frame(leverage = leverage,
                                     residuals = residuals_vals),
                   aes(x = leverage, y = residuals)) +
        geom_point(alpha = 0.6, size = 2) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        labs(title = paste("Residuals vs Leverage -", response_var),
             x = "Leverage", y = "Residuals") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
    } else {
      p4 <- ggplot() +
        theme_void() +
        labs(title = "Leverage plot not available for mixed models")
    }
    
    plots <- list(
      residuals_fitted = p1,
      qq_plot = p2,
      scale_location = p3,
      leverage = p4
    )
    
    #Statistical tests (for lm)
    if (inherits(model, "lm")) {
      #Normality test
      normality_test <- shapiro.test(residuals_vals)
      diagnostics$normality_test <- normality_test
      
      #Heteroscedasticity test
      hetero_test <- car::ncvTest(model)
      diagnostics$heteroscedasticity_test <- hetero_test
      
      cat("Model diagnostics for", response_var, ":\n")
      cat("  - Normality (Shapiro-Wilk): W =", round(normality_test$statistic, 4),
          "p =", round(normality_test$p.value, 4), "\n")
      cat("  - Heteroscedasticity (NCV): Chi² =", round(hetero_test$ChiSquare, 4),
          "p =", round(hetero_test$p, 4), "\n")
    }
  }
  
  return(list(plots = plots, diagnostics = diagnostics))
}

#ENHANCED CORRELATION MATRIX FUNCTION
# ====================================
create_correlation_plot <- function(data, response_vars, method = "pearson") {
  #Select only numeric variables for correlation
  numeric_data <- data %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(any_of(response_vars)) %>%
    na.omit()
  
  if (ncol(numeric_data) < 2) {
    cat("Insufficient numeric variables for correlation matrix\n")
    return(NULL)
  }
  
  cor_matrix <- cor(numeric_data, method = method)
  
  p <- ggcorrplot::ggcorrplot(
    cor_matrix,
    method = "circle",
    type = "lower",
    lab = TRUE,
    lab_size = 3,
    colors = c("#6D9EC1", "white", "#E46726"),
    outline.color = "white",
    show.legend = TRUE,
    title = "Correlation Matrix of Response Variables",
    ggtheme = theme_minimal()
  ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

#NON-PARAMETRIC ANALYSIS FUNCTION
#Kruskal–Wallis across groups and BH-adjusted pairwise Wilcoxon.
# 'exact = FALSE' avoids tie warnings; results remain valid (large-sample approx).
nonparametric_analysis <- function(data, response_vars, group_var = "Location") {
  results <- list()
  
  cat("\n=== NON-PARAMETRIC ANALYSIS ===\n")
  
  for (resp in response_vars) {
    cat("\n--- Analysis for:", resp, "---\n")
    
    # Kruskal-Wallis test for overall group differences
    kruskal_formula <- as.formula(paste(resp, "~", group_var))
    kruskal_test <- kruskal.test(kruskal_formula, data = data)
    
    cat("Kruskal-Wallis Test:\n")
    cat("  H =", round(kruskal_test$statistic, 3),
        "df =", kruskal_test$parameter,
        "p =", round(kruskal_test$p.value, 4), "\n")
    
    #Pairwise Wilcoxon tests with p-value adjustment (ties-safe)
    if (kruskal_test$p.value < 0.05) {
      pairwise_wilcox <- pairwise.wilcox.test(
        data[[resp]], data[[group_var]],
        p.adjust.method = "BH",
        exact = FALSE,            # <-- avoids exact-p with ties
        na.action = na.omit
      )
      cat("Pairwise Wilcoxon Tests (BH-adjusted):\n")
      print(pairwise_wilcox$p.value)
    } else {
      pairwise_wilcox <- NULL
      cat("No significant overall effect - skipping pairwise tests\n")
    }
    
    #Effect size - epsilon squared
    n <- nrow(data)
    k <- length(unique(data[[group_var]]))
    epsilon_squared <- as.numeric((kruskal_test$statistic - (k - 1)) / (n - k))
    
    results[[resp]] <- list(
      kruskal_wallis = kruskal_test,
      pairwise_wilcoxon = pairwise_wilcox,
      effect_size = epsilon_squared,
      n_groups = k,
      total_n = n
    )
  }
  
  return(results)
}


#COMPREHENSIVE PCA PLOTTING FUNCTION
#  Performs PCA on numeric response matrix, returns (1) scores plot with
# 95% ellipses (only for groups with n >= 3), (2) variance explained, and
# (3) loadings plot. Ellipse condition prevents "no non-missing arguments" warnings.
create_pca_plots_comprehensive <- function(pca_result, metadata, response_vars) {
  pca_scores <- as.data.frame(pca_result$x)
  if (!is.null(metadata$Location)) {
    pca_scores$Location <- metadata$Location
  } else if ("Location" %in% colnames(metadata)) {
    pca_scores$Location <- metadata$Location
  } else {
    pca_scores$Location <- factor("Unknown")
  }
  variance_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
  
  #Only draw ellipses for groups with >= 3 points
  group_sizes <- table(pca_scores$Location)
  has_ellipse <- names(group_sizes[group_sizes >= 3])
  
  #1.Main PCA plot with conditional ellipses
  pca_main <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Location)) +
    geom_point(size = 3, alpha = 0.8) +
    { if (length(has_ellipse) > 0)
      stat_ellipse(data = subset(pca_scores, Location %in% has_ellipse),
                   level = 0.95, alpha = 0.2, size = 1) else NULL } +
    labs(
      x = paste0("PC1 (", variance_explained[1], "%)"),
      y = paste0("PC2 (", variance_explained[2], "%)"),
      title = "Principal Component Analysis",
      subtitle = "Colored by Location (ellipses shown only for groups with n ≥ 3)",
      color = "Location"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_color_viridis_d(option = "plasma")
  
  #2.Variance explained plot
  variance_data <- data.frame(
    PC = paste0("PC", 1:length(variance_explained)),
    Variance = variance_explained,
    Cumulative = cumsum(variance_explained)
  )
  variance_data <- variance_data[1:min(8, nrow(variance_data)), ]
  
  pca_variance <- ggplot(variance_data, aes(x = reorder(PC, -Variance), y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = Cumulative, group = 1), color = "red", size = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 2) +
    labs(
      title = "PCA Variance Explained",
      subtitle = "Bar: Individual PC, Line: Cumulative",
      x = "Principal Component",
      y = "Variance Explained (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(sec.axis = sec_axis(~./100, name = "Cumulative Proportion"))
  
  #3.Loadings plot
  loadings <- as.data.frame(pca_result$rotation[, 1:2, drop = FALSE])
  loadings$Variable <- rownames(loadings)
  
  pca_loadings <- ggplot(loadings, aes(x = PC1, y = PC2)) +
    geom_segment(aes(xend = 0, yend = 0),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "darkred", alpha = 0.7) +
    geom_text(aes(label = Variable),
              hjust = -0.1, vjust = -0.1,
              size = 3, fontface = "bold") +
    labs(
      title = "PCA Variable Loadings",
      subtitle = "Arrows show variable contributions to PCs",
      x = paste0("PC1 (", variance_explained[1], "%)"),
      y = paste0("PC2 (", variance_explained[2], "%)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  return(list(
    main = pca_main,
    variance = pca_variance,
    loadings = pca_loadings
  ))
}


#ANALYSIS CONFIGURATION SYSTEM
#########
###
##pare5BROF,pare5BRIF,pare5LR,pare5BRin,pare5BRout,
#resid15BROF,resid15BRIF,resid15LR,resid15BRin,resid15BRout,
#transi15BROF,transi15BRIF,transi15LR,transi15BRin,transi15BRout
####
#Clean_Data$logLR<-log1p(Clean_Data$crypto5LR)
#Clean_Data$logIF<-log1p(Clean_Data$resid15BRIF)
#Clean_Data$netflow<- (Clean_Data$crypto5BROF - Clean_Data$crypto5BRIF)/(Clean_Data$crypto5BROF + Clean_Data$crypto5BRIF)




#########
analysis_config <- list(
  core_ecological = list(
    response_vars = c("biomass_kg_ha", "abundance_ind_250m2", "FRic", "FEve"),
    predictor_vars = c("Geomorphology", "sedimnt", "nutrint","NO_TK_AREA" ,"NO_TK_AREA" , "Location","crypto5BROF","crypto5BRIF","crypto5LR","crypto5BRin","crypto5BRout","sst_sd_6_year_mean","sst_q90_6_year_mean","sst_mean_6_year_mean"),
    random_effects = c("(1|Location)"),
    description = "Core ecological variables analysis"
  ),
  functional_diversity = list(
    response_vars = c("FRic", "FEve", "FDiv", "FDis", "RaoQ"),
    predictor_vars = c("Geomorphology", "sedimnt", "nutrint","NO_TK_AREA" ,"NO_TK_AREA" , "Location","crypto5BROF","crypto5BRIF","crypto5LR","crypto5BRin","crypto5BRout","sst_sd_6_year_mean","sst_q90_6_year_mean","sst_mean_6_year_mean"),
    random_effects = c("(1|Location)"),
    description = "Functional diversity metrics analysis"
  ),
  biomass_focused = list(
    response_vars = c("biomass_kg_ha", "abundance_ind_250m2"),
    predictor_vars = c("Geomorphology", "sedimnt", "nutrint","NO_TK_AREA" ,"NO_TK_AREA" , "Location","crypto5BROF","crypto5BRIF","crypto5LR","crypto5BRin","crypto5BRout","sst_sd_6_year_mean","sst_q90_6_year_mean","sst_mean_6_year_mean"),
    random_effects = c("(1|Location)"),
    description = "Biomass and abundance focused analysis"
  )
)



#ANALYSIS FUNCTION
#  descriptive stats, plots, non-parametrics,
# linear models, mixed models (with optional Location-as-random-only toggle),
# PCA, PERMANOVA, and saving of all artifacts.
analyze_reef_data_enhanced <- function(
    response_vars = c("biomass_kg_ha", "abundance_ind_250m2", "FRic"),
    predictor_vars = c("Geomorphology", "sedimnt", "nutrint", "Location"),
    random_effects = c("(1|Location)"),
    data = clean_data,
    output_prefix = "enhanced_analysis",
    create_plots = TRUE,
    save_results = TRUE,
    run_diagnostics = TRUE,
    run_nonparametric = TRUE
) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("STARTING ENHANCED REEF ECOLOGY ANALYSIS\n")
  cat(strrep("=", 70), "\n")
  
  # --- NEW: toggle to avoid Location as both fixed and random (keeps default as-is)
  lmm_location_as_random_only <- TRUE  # set TRUE to drop 'Location' from fixed effects when (1|Location) used
  
  # Input validation
  if (length(response_vars) == 0) stop("No response variables specified.")
  if (length(predictor_vars) == 0) stop("No predictor variables specified.")
  
  # Enhanced data validation
  validation_issues <- validate_analysis_data(data, response_vars, predictor_vars)
  if (length(validation_issues$errors) > 0) {
    stop("Critical data validation errors detected. Analysis aborted.")
  }
  
  # Create analysis dataset (any_of avoids hard errors if a column is absent)
  analysis_data <- data %>%
    dplyr::select(dplyr::any_of(c(
      response_vars, predictor_vars, "Site_ID", "Transect", "Dive_Site", "Location"
    ))) %>%
    dplyr::filter(stats::complete.cases(.))
  
  if (nrow(analysis_data) == 0) {
    stop("No complete cases after filtering. Check variable names and missing data.")
  }
  
  #light scaling for continuous predictors to improve conditioning
  #(does not alter saved raw summaries; only improves model fit numerics)
  num_cols <- intersect(c("sedimnt", "nutrint"), names(analysis_data))
  if (length(num_cols) > 0) {
    analysis_data[num_cols] <- lapply(analysis_data[num_cols], function(z) as.numeric(scale(z)))
  }
  
  cat("Analysis dataset created:\n")
  cat("  - Variables:", length(unique(c(response_vars, predictor_vars))), "total\n")
  cat("  - Observations:", nrow(analysis_data), "complete cases\n")
  locs <- unique(analysis_data$Location)
  cat("  - Locations:", paste(utils::head(locs, 20), collapse = ", "),
      if (length(locs) > 20) " ... [truncated]" else "", "\n", sep = "")
  
  results <- list(
    metadata = list(
      response_vars   = response_vars,
      predictor_vars  = predictor_vars,
      random_effects  = random_effects,
      sample_size     = nrow(analysis_data),
      locations       = locs,
      n_locations     = length(locs),
      timestamp       = Sys.time(),
      validation_issues = validation_issues
    )
  )
  
  #1.DESCRIPTIVE STATISTICS
  # ==================================
  cat("\n=== ENHANCED DESCRIPTIVE STATISTICS ===\n")
  
  location_summary <- analysis_data %>%
    dplyr::group_by(.data$Location) %>%
    dplyr::summarise(
      n_transects = dplyr::n(),
      dplyr::across(
        dplyr::all_of(response_vars),
        list(
          mean   = ~ mean(.x, na.rm = TRUE),
          sd     = ~ stats::sd(.x, na.rm = TRUE),
          median = ~ stats::median(.x, na.rm = TRUE),
          min    = ~ min(.x, na.rm = TRUE),
          max    = ~ max(.x, na.rm = TRUE),
          cv     = ~ ifelse(mean(.x, na.rm = TRUE) == 0, NA_real_,
                            stats::sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE))
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  
  overall_summary <- analysis_data %>%
    dplyr::summarise(
      n_total     = dplyr::n(),
      n_locations = dplyr::n_distinct(.data$Location),
      dplyr::across(
        dplyr::all_of(response_vars),
        list(
          mean   = ~ mean(.x, na.rm = TRUE),
          sd     = ~ stats::sd(.x, na.rm = TRUE),
          median = ~ stats::median(.x, na.rm = TRUE),
          min    = ~ min(.x, na.rm = TRUE),
          max    = ~ max(.x, na.rm = TRUE),
          cv     = ~ ifelse(mean(.x, na.rm = TRUE) == 0, NA_real_,
                            stats::sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE))
        ),
        .names = "{.col}_{.fn}"
      )
    )
  
  print(location_summary)
  
  results[["descriptive"]] <- list(
    by_location = location_summary,
    overall     = overall_summary
  )
  
  if (isTRUE(save_results)) {
    dir.create("tables", showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(location_summary,
                     file = file.path("tables", paste0(output_prefix, "_location_summary.csv")),
                     row.names = FALSE)
    utils::write.csv(overall_summary,
                     file = file.path("tables", paste0(output_prefix, "_overall_summary.csv")),
                     row.names = FALSE)
  }
  
  #2. VISUALIZATION

  plots <- list()
  if (isTRUE(create_plots)) {
    cat("\n=== ENHANCED DATA VISUALIZATION ===\n")
    
    #2.1 Boxplots for each response variable
    for (resp in response_vars) {
      if (!resp %in% names(analysis_data)) next
      p_box <- ggplot(analysis_data, aes(x = Location, y = .data[[resp]], fill = Location)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
        stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
        labs(title = paste("Distribution of", resp, "by Location"),
             subtitle = paste("n =", nrow(analysis_data), "transects; Red diamond = mean"),
             y = resp,
             x = "Location") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none",
              plot.title = element_text(face = "bold")) +
        scale_fill_viridis_d(option = "viridis")
      
      plots[[paste("boxplot", resp, sep = "_")]] <- p_box
      print(p_box)
      
      if (isTRUE(save_results)) {
        ggsave(file.path("figures", paste0(output_prefix, "_boxplot_", resp, ".png")),
               p_box, width = 10, height = 6, dpi = 300, bg = "white")
      }
    }
    
    #2.2 Correlation matrix
    cor_plot <- create_correlation_plot(analysis_data, response_vars)
    if (!is.null(cor_plot)) {
      plots[["correlation_matrix"]] <- cor_plot
      print(cor_plot)
      if (isTRUE(save_results)) {
        ggsave(file.path("figures", paste0(output_prefix, "_correlation_matrix.png")),
               cor_plot, width = 8, height = 6, dpi = 300, bg = "white")
      }
    }
  }
  results[["plots"]] <- plots
  
  #3.NON-PARAMETRIC ANALYSIS

  if (isTRUE(run_nonparametric)) {
    nonparametric_results <- nonparametric_analysis(analysis_data, response_vars)
    results[["nonparametric"]] <- nonparametric_results
    
    if (isTRUE(save_results)) {
      nonparametric_df <- purrr::map_dfr(names(nonparametric_results), function(resp) {
        result <- nonparametric_results[[resp]]
        data.frame(
          Response = resp,
          Test = "Kruskal-Wallis",
          Statistic = as.numeric(result$kruskal_wallis$statistic),
          df = as.numeric(result$kruskal_wallis$parameter),
          p_value = as.numeric(result$kruskal_wallis$p.value),
          Effect_Size = as.numeric(result$effect_size),
          n_groups = result$n_groups,
          total_n = result$total_n
        )
      })
      write.csv(nonparametric_df,
                file.path("tables", paste0(output_prefix, "_nonparametric_results.csv")),
                row.names = FALSE)
    }
  }
  
  #4.LINEAR MODELS WITH ENHANCED DIAGNOSTICS

  cat("\n=== LINEAR MODELS WITH ENHANCED DIAGNOSTICS ===\n")
  
  lm_results <- list()
  lm_diagnostics <- list()
  
  for (resp in response_vars) {
    if (!resp %in% names(analysis_data)) next
    cat("\n--- Linear Model for:", resp, "---\n")
    
    lm_formula <- as.formula(paste(resp, "~", paste(predictor_vars, collapse = " + ")))
    lm_model <- lm(lm_formula, data = analysis_data)
    
    lm_summary <- summary(lm_model)
    cat("Model Summary:\n"); print(lm_summary)
    
    anova_result <- tryCatch({
      car::Anova(lm_model, type = "III")
    }, error = function(e) {
      cat("Type III ANOVA failed, trying Type II:\n")
      tryCatch({
        car::Anova(lm_model, type = "II")
      }, error = function(e2) {
        cat("All ANOVA types failed, using basic anova\n")
        anova(lm_model)
      })
    })
    
    if (!is.null(anova_result)) {
      cat("\nANOVA Results:\n"); print(anova_result)
    }
    
    if (isTRUE(run_diagnostics)) {
      model_diagnostics <- create_model_diagnostics(lm_model, resp, "Linear Model")
      lm_diagnostics[[resp]] <- model_diagnostics
      
      if (isTRUE(save_results) && isTRUE(create_plots)) {
        diagnostic_plots <- model_diagnostics$plots
        for (plot_name in names(diagnostic_plots)) {
          ggsave(file.path("diagnostics", paste0(output_prefix, "_lm_", resp, "_", plot_name, ".png")),
                 diagnostic_plots[[plot_name]], width = 8, height = 6, dpi = 300, bg = "white")
        }
      }
    }
    
    #Multicollinearity screening (VIF) for LM + inline flagged summary
    #Uses car::vif; for multi-df factors, we report GVIF and GVIF^(1/(2*Df)).
    vif_tbl <- NULL
    vif_tbl <- tryCatch({
      vf <- car::vif(lm_model)
      if (is.matrix(vf)) {
        out <- as.data.frame(vf)
        out$Term <- rownames(out)
        out$GVIF_adj <- out$GVIF^(1/(2*out$Df))
        rownames(out) <- NULL
        out <- out[, c("Term","Df","GVIF","GVIF_adj")]
      } else {
        out <- data.frame(
          Term = names(vf),
          Df = 1L,
          GVIF = as.numeric(vf),
          GVIF_adj = as.numeric(vf) # same for 1 df
        )
      }
    }, error = function(e) {
      message("VIF computation failed for ", resp, ": ", conditionMessage(e))
      NULL
    })
    
    # Inline summary + threshold flags
    if (!is.null(vif_tbl)) {
      # Flags: GVIF > 5 (classic) OR GVIF_adj > 2 (rough heuristic for multi-df terms)
      vif_tbl <- vif_tbl %>%
        mutate(Flag = (GVIF > 5) | (GVIF_adj > 2))
      n_flags <- sum(vif_tbl$Flag, na.rm = TRUE)
      cat("\nVIF Summary (", resp, "):\n", sep = "")
      print(vif_tbl)
      if (n_flags > 0) {
        cat("⚠️  Multicollinearity flags detected in", n_flags, "term(s):\n")
        print(vif_tbl %>% dplyr::filter(Flag) %>% dplyr::select(Term, Df, GVIF, GVIF_adj))
        cat("Note: Consider interpretation caution, reparameterization, or removing/re-grouping highly collinear terms.\n")
      } else {
        cat("No VIF flags (GVIF > 5 or GVIF_adj > 2) detected.\n")
      }
      if (isTRUE(save_results)) {
        dir.create("tables", showWarnings = FALSE, recursive = TRUE)
        utils::write.csv(vif_tbl,
                         file = file.path("tables", paste0(output_prefix, "_VIF_", resp, ".csv")),
                         row.names = FALSE)
      }
    }
    
    lm_results[[resp]] <- list(
      model = lm_model,
      summary = lm_summary,
      anova = anova_result,
      coefficients = broom::tidy(lm_model),
      performance = broom::glance(lm_model),
      diagnostics = if (isTRUE(run_diagnostics)) lm_diagnostics[[resp]] else NULL,
      vif = vif_tbl # --- NEW: persist VIF table per response
    )
  }
  
  results[["linear_models"]] <- lm_results
  
  if (isTRUE(save_results) && length(lm_results) > 0) {
    lm_coefs <- purrr::map_dfr(lm_results, function(x) x$coefficients, .id = "Response")
    write.csv(lm_coefs,
              file.path("tables", paste0(output_prefix, "_linear_model_coefficients.csv")),
              row.names = FALSE)
    lm_performance <- purrr::map_dfr(lm_results, function(x) x$performance, .id = "Response")
    write.csv(lm_performance,
              file.path("tables", paste0(output_prefix, "_linear_model_performance.csv")),
              row.names = FALSE)
    
    #Combined VIF across responses
    all_vif <- purrr::map_dfr(names(lm_results), function(nm) {
      vt <- lm_results[[nm]]$vif
      if (is.null(vt)) return(NULL)
      cbind(Response = nm, vt)
    })
    if (nrow(all_vif) > 0) {
      utils::write.csv(
        all_vif,
        file = file.path("tables", paste0(output_prefix, "_VIF_all_responses.csv")),
        row.names = FALSE
      )
    }
  }
  
  #5.MIXED EFFECTS MODELS WITH checks

  cat("\n=== MIXED EFFECTS MODELS WITH ENHANCED DIAGNOSTICS ===\n")
  
  mixed_results <- list()
  mixed_diagnostics <- list()
  
  for (resp in response_vars) {
    if (!resp %in% names(analysis_data)) next
    cat("\n--- Mixed Model for:", resp, "---\n")
    
    # Optionally remove 'Location' from fixed effects if also used as random intercept
    fixed_terms <- paste(predictor_vars, collapse = " + ")
    if (lmm_location_as_random_only) {
      fixed_terms <- paste(setdiff(predictor_vars, "Location"), collapse = " + ")
      fixed_terms <- if (fixed_terms == "") "1" else fixed_terms
    }
    
    mixed_formula <- as.formula(paste(
      resp, "~", fixed_terms, "+", paste(random_effects, collapse = " + ")
    ))
    
    mixed_model <- tryCatch({
      lme4::lmer(mixed_formula, data = analysis_data,
                 control = lme4::lmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 100000)))
    }, error = function(e) {
      cat("Complex mixed model failed, trying simpler structure:\n")
      simple_formula <- as.formula(paste(
        resp, "~", fixed_terms, "+ (1|Location)"
      ))
      lme4::lmer(simple_formula, data = analysis_data,
                 control = lme4::lmerControl(optimizer = "bobyqa"))
    })
    
    if (!inherits(mixed_model, "try-error") && !is.null(mixed_model)) {
      mixed_summary <- summary(mixed_model)
      cat("Mixed Model Summary:\n"); print(mixed_summary)
      
      r2 <- tryCatch({
        MuMIn::r.squaredGLMM(mixed_model)
      }, error = function(e) {
        cat("R-squared calculation failed\n"); NULL
      })
      if (!is.null(r2)) {
        cat("\nR-squared (Marginal + Conditional):\n"); print(r2)
      }
      
      variance_components <- as.data.frame(VarCorr(mixed_model))
      cat("\nVariance Components:\n"); print(variance_components)
      
      if (isTRUE(run_diagnostics)) {
        dh <- tryCatch({
          DHARMa::simulateResiduals(mixed_model, n = 1000)
        }, error = function(e) NULL)
        if (!is.null(dh)) {
          mixed_diagnostics[[resp]] <- list(dharma = dh)
          if (isTRUE(save_results) && isTRUE(create_plots)) {
            png(file.path("diagnostics", paste0(output_prefix, "_lmer_", resp, "_DHARMa.png")),
                width = 900, height = 900, res = 150)
            plot(dh)
            dev.off()
          }
        }
      }
      
      mixed_results[[resp]] <- list(
        model = mixed_model,
        summary = mixed_summary,
        r2 = r2,
        varcomp = variance_components,
        diagnostics = if (isTRUE(run_diagnostics)) mixed_diagnostics[[resp]] else NULL
      )
    } else {
      cat("Mixed model failed for response:", resp, "\n")
    }
  }
  
  results[["mixed_models"]] <- mixed_results
  
  if (isTRUE(save_results) && length(mixed_results) > 0) {
    sink(file.path("tables", paste0(output_prefix, "_mixed_model_summaries.txt")))
    for (nm in names(mixed_results)) {
      cat("\n====== Mixed Model:", nm, "======\n")
      print(mixed_results[[nm]]$summary)
      if (!is.null(mixed_results[[nm]]$r2)) {
        cat("\nR2 (marginal/conditional):\n"); print(mixed_results[[nm]]$r2)
      }
      cat("\nVariance Components:\n"); print(mixed_results[[nm]]$varcomp)
    }
    sink()
    lapply(names(mixed_results), function(nm) {
      saveRDS(mixed_results[[nm]]$model, file.path("models", paste0(output_prefix, "_lmer_", nm, ".rds")))
    })
  }
  
  # --- NEW: MODEL COMPARISON (AIC/BIC) BETWEEN LM AND LMM --------------------
  # Rationale: compare information criteria under the same likelihood basis.
  # We keep LMMs fitted via REML for inference; refit with ML (REML=FALSE) for IC comparison.
  cat("\n=== MODEL COMPARISON (AIC/BIC): LM vs LMM ===\n")
  model_comp_rows <- list()
  if (length(lm_results) > 0 && length(mixed_results) > 0) {
    for (resp in names(lm_results)) {
      lm_obj <- lm_results[[resp]]$model
      lmm_obj <- if (!is.null(mixed_results[[resp]])) mixed_results[[resp]]$model else NULL
      if (is.null(lm_obj) || is.null(lmm_obj)) {
        cat("No comparable pair (LM/LMM) for:", resp, "\n")
        next
      }
      # Refit LMM with ML for fair comparison
      lmm_ml <- tryCatch({
        update(lmm_obj, REML = FALSE)
      }, error = function(e) {
        cat("ML refit failed for response:", resp, " | ", conditionMessage(e), "\n")
        NULL
      })
      if (is.null(lmm_ml)) next
      
      aic_lm  <- AIC(lm_obj);  bic_lm  <- BIC(lm_obj)
      aic_lmm <- AIC(lmm_ml);  bic_lmm <- BIC(lmm_ml)
      delta_aic <- aic_lmm - aic_lm
      delta_bic <- bic_lmm - bic_lm
      
      cat(sprintf("  %s: LM AIC=%.2f, LMM(ML) AIC=%.2f, ΔAIC=%.2f | LM BIC=%.2f, LMM(ML) BIC=%.2f, ΔBIC=%.2f\n",
                  resp, aic_lm, aic_lmm, delta_aic, bic_lm, bic_lmm, delta_bic))
      
      model_comp_rows[[resp]] <- data.frame(
        Response = resp,
        LM_AIC = aic_lm,
        LMM_AIC = aic_lmm,
        Delta_AIC = delta_aic,  # negative -> LMM better
        LM_BIC = bic_lm,
        LMM_BIC = bic_lmm,
        Delta_BIC = delta_bic   # negative -> LMM better
      )
    }
  }
  model_comp_df <- dplyr::bind_rows(model_comp_rows)
  results[["model_comparison"]] <- model_comp_df
  if (isTRUE(save_results) && nrow(model_comp_df) > 0) {
    utils::write.csv(
      model_comp_df,
      file = file.path("tables", paste0(output_prefix, "_model_comparison_AIC_BIC.csv")),
      row.names = FALSE
    )
  }
  # --------------------------------------------------------------------------
  
  # 6. PCA & PERMANOVA
  # ===================
  cat("\n=== MULTIVARIATE ANALYSIS: PCA & PERMANOVA ===\n")
  
  resp_mat <- analysis_data %>%
    dplyr::select(any_of(response_vars)) %>%
    dplyr::select(where(is.numeric)) %>%
    as.data.frame()
  
  if (ncol(resp_mat) >= 2) {
    pca_res <- prcomp(resp_mat, center = TRUE, scale. = TRUE)
    pca_plots <- create_pca_plots_comprehensive(pca_res, metadata = analysis_data, response_vars = response_vars)
    results[["pca"]] <- list(result = pca_res, plots = pca_plots)
    
    if (isTRUE(create_plots)) {
      print(pca_plots$main); print(pca_plots$variance); print(pca_plots$loadings)
    }
    if (isTRUE(save_results)) {
      ggsave(file.path("ordination", paste0(output_prefix, "_PCA_main.png")),
             pca_plots$main, width = 9, height = 7, dpi = 300, bg = "white")
      ggsave(file.path("ordination", paste0(output_prefix, "_PCA_variance.png")),
             pca_plots$variance, width = 9, height = 6, dpi = 300, bg = "white")
      ggsave(file.path("ordination", paste0(output_prefix, "_PCA_loadings.png")),
             pca_plots$loadings, width = 9, height = 7, dpi = 300, bg = "white")
      saveRDS(pca_res, file.path("ordination", paste0(output_prefix, "_PCA_result.rds")))
    }
  } else {
    cat("PCA skipped: need at least two numeric response variables.\n")
  }
  
  permanova_overall <- NULL
  pairwise_perm <- NULL
  if (ncol(resp_mat) >= 1 && "Location" %in% names(analysis_data)) {
    dist_matrix <- vegan::vegdist(resp_mat, method = "euclidean")
    permanova_overall <- vegan::adonis2(dist_matrix ~ Location, data = analysis_data, permutations = 999)
    cat("\nOverall PERMANOVA by Location:\n"); print(permanova_overall)
    
    pairwise_perm <- pairwise.adonis2(resp_mat, analysis_data, permutations = 999, method = "euclidean", p_adjust = "BH")
    cat("\nPairwise PERMANOVA (BH-adjusted):\n"); print(pairwise_perm$pairwise_results)
    
    results[["permanova"]] <- list(overall = permanova_overall, pairwise = pairwise_perm)
    if (isTRUE(save_results)) {
      sink(file.path("tables", paste0(output_prefix, "_permanova_overall.txt"))); print(permanova_overall); sink()
      if (nrow(pairwise_perm$pairwise_results) > 0) {
        write.csv(pairwise_perm$pairwise_results,
                  file.path("tables", paste0(output_prefix, "_permanova_pairwise.csv")),
                  row.names = FALSE)
      }
    }
  } else {
    cat("PERMANOVA skipped: insufficient response matrix or missing Location.\n")
  }
  
  # 7. SAVE MODELS (LM) AS RDS
  # ===========================
  if (isTRUE(save_results) && length(lm_results) > 0) {
    lapply(names(lm_results), function(nm) {
      saveRDS(lm_results[[nm]]$model, file.path("models", paste0(output_prefix, "_lm_", nm, ".rds")))
    })
  }
  
  # 8. FINAL REPORT STUB
  # =====================
  report_stub <- list(
    note = "Report generation placeholder. Integrate with R Markdown if needed.",
    created = Sys.time(),
    outputs = list(
      figures = list.files("figures", full.names = TRUE),
      tables  = list.files("tables", full.names = TRUE),
      models  = list.files("models", full.names = TRUE),
      ordination = list.files("ordination", full.names = TRUE),
      diagnostics = list.files("diagnostics", full.names = TRUE)
    )
  )
  results[["report"]] <- report_stub
  
  cat("\n", strrep("=", 70), "\n")
  cat("ANALYSIS COMPLETE. RESULTS OBJECT RETURNED.\n")
  cat(strrep("=", 70), "\n")
  
  return(results)
}


# =========================
# EXAMPLE RUN(S) (EDIT ME)
# =========================
# You can run any of the predefined configs or provide your own variables.
# Uncomment ONE of the following if you want the script to run on source.

cfg <- analysis_config$core_ecological
res <- analyze_reef_data_enhanced(
  response_vars = cfg$response_vars,
  predictor_vars = cfg$predictor_vars,
  random_effects = cfg$random_effects,
  data = clean_data,
  output_prefix = "enhanced_core_ecological",
  create_plots = TRUE,
  save_results = TRUE,
  run_diagnostics = TRUE,
  run_nonparametric = TRUE
)

cfg <- analysis_config$functional_diversity
res_fd <- analyze_reef_data_enhanced(
  response_vars = cfg$response_vars,
  predictor_vars = cfg$predictor_vars,
  random_effects = cfg$random_effects,
  data = clean_data,
  output_prefix = "enhanced_functional_diversity",
  create_plots = TRUE,
  save_results = TRUE,
  run_diagnostics = TRUE,
  run_nonparametric = TRUE
)

cfg <- analysis_config$biomass_focused
res_biomass <- analyze_reef_data_enhanced(
  response_vars = cfg$response_vars,
  predictor_vars = cfg$predictor_vars,
  random_effects = cfg$random_effects,
  data = clean_data,
  output_prefix = "enhanced_biomass_focused",
  create_plots = TRUE,
  save_results = TRUE,
  run_diagnostics = TRUE,
  run_nonparametric = TRUE
)

