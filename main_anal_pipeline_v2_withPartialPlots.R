
#ANALYSIS PIPELINE  — FIji Analyses
#Data validation & cleaning
#Linear (LM) & Mixed (LMM) models with collinearity control
#Predictor scaling (LMM always; LM configurable)
#Univariate panel for excluded predictors (LM + LMM)
#Non-parametric tests (Kruskal–Wallis, pairwise Wilcoxon BH)
#PCA & PERMANOVA
#Diagnostics, figures, tables, RDS models, report stub
#Type III tests for LMMs; singular-LMM auto-flag; corrected VIF logic
#Partial plots for LM and LMM
#======================================================================
rm(list = ls())
cat("Initializing Reef Ecology Analysis Pipeline...\n")

# ---------------------------
#Package install and loading
# ---------------------------
install_and_load_packages <- function() {
  required_packages <- c(
    "tidyverse","lme4","lmerTest","vegan","ggplot2",
    "car","MuMIn","DHARMa","ggpubr","GGally",
    "broom","broom.mixed","glue","ggcorrplot","patchwork","sandwich","lmtest", "clubSandwich","pbkrtest"
  )
  inst <- rownames(installed.packages())
  missing_packages <- setdiff(required_packages, inst)
  if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, dependencies = TRUE)
  }
  cat("Loading required packages...\n")
  suppressPackageStartupMessages({
    library(tidyverse); library(lme4); library(lmerTest); library(vegan)
    library(ggplot2);  library(car);  library(MuMIn);   library(DHARMa)
    library(ggpubr);   library(GGally); library(broom); library(broom.mixed)
    library(glue);     library(ggcorrplot); library(patchwork);library(sandwich)
    library(lmtest);library(clubSandwich); library(pbkrtest)
  })
  cat("ll packages loaded successfully!\n")
}
install_and_load_packages()

#Sumtozero contrasts for valid Type III tests
options(contrasts = c("contr.sum", "contr.poly"))

#onetime init for model chooser table
choice_tbl <- data.frame(
  response = character(),
  AIC_LM = double(),
  AIC_LMM = double(),
  delta = double(),
  singular = logical(),
  choice = character(),
  stringsAsFactors = FALSE
)

#Helper forType III with a specific df method
attempt_type3 <- function(model, ddf = c("Kenward-Roger","Satterthwaite")) {
  ddf <- match.arg(ddf)
  tryCatch(
    anova(model, type = 3, ddf = ddf),  # rely on S3 dispatch from lmerTest
    error = function(e) { message("Type III (", ddf, ") failed: ", conditionMessage(e)); NULL }
  )
}


#For reproducible steps (screening, DHARMa sims etc)
set.seed(42)

#Working dir and  data loading

cat("Setting up working directory and loading data...\n")
if (!exists("custom_wd")) {
  custom_wd <- getwd()
  cat("Using working directory:", custom_wd, "\n")
}
if (!dir.exists(custom_wd)) stop("Working directory does not exist: ", custom_wd)
setwd(custom_wd)

load_reef_data <- function(file_path = "DOV_transect_connect_final.csv") {
  if (!file.exists(file_path)) {
    stop("Data file not found: ", file_path,
         "\nPlease ensure the file exists in the working directory.")
  }
  cat("Loading data from:", file_path, "\n")
  data <- read.csv(file_path)
  cat("Data loaded successfully:", nrow(data), "rows,", ncol(data), "columns\n")
  data
}
data <- load_reef_data()

#explore the dataset to see what vaeiables have bene loaded and think anout your reapose variables
# Output directiries
create_output_directories <- function() {
  dirs <- c("figures","tables","models","ordination","diagnostics","reports")
  for (dir in dirs) if (!dir.exists(dir)) dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  cat("output directories ready.\n")
}
create_output_directories()

purge_outputs_by_prefix <- function(prefix) {
  pats <- c("figures","tables","models","ordination","diagnostics")
  for (d in pats) {
    files <- list.files(d, pattern = paste0("^", prefix, "_"), full.names = TRUE)
    if (length(files)) file.remove(files)
  }
}

#check the data for formats etc

validate_analysis_data <- function(data, response_vars, predictor_vars) {
  cat("=== DATA VALIDATION CHECKS ===\n")
  issues <- list(warnings = character(), errors = character())
  if (nrow(data) == 0) {
    issues$errors <- c(issues$errors, "Dataset is empty")
    return(issues)
  }
  missing_vars <- c(setdiff(response_vars, names(data)),
                    setdiff(predictor_vars, names(data)))
  if (length(missing_vars) > 0) {
    issues$errors <- c(issues$errors, paste("Missing variables:", paste(missing_vars, collapse = ", ")))
  }
  if (nrow(data) < 10) {
    issues$warnings <- c(issues$warnings, "Small sample size (n < 10) may affect model reliability")
  }
  for (var in response_vars) if (var %in% names(data)) {
    vv <- data[[var]]
    if (all(is.na(vv))) issues$warnings <- c(issues$warnings, paste("All NA in response:", var))
    if (var(vv, na.rm = TRUE) == 0) issues$warnings <- c(issues$warnings, paste("Zero variance response:", var))
    denom <- sum(!is.na(vv)); if (denom > 0) {
      zero_prop <- sum(vv == 0, na.rm = TRUE) / denom
      if (zero_prop > 0.8) issues$warnings <- c(issues$warnings, glue("High zeros ({round(zero_prop*100,1)}%) in: {var}"))
    }
  }
  existing_preds <- predictor_vars[predictor_vars %in% names(data)]
  facs <- existing_preds[sapply(existing_preds, \(x) is.factor(data[[x]]))]
  for (f in facs) {
    if (length(unique(na.omit(data[[f]]))) < 2) issues$errors <- c(issues$errors, paste("Factor has 1 level:", f))
  }
  cc_cols <- c(intersect(response_vars, names(data)), intersect(predictor_vars, names(data)))
  if (length(cc_cols) > 0) {
    complete_cases <- complete.cases(data[, cc_cols, drop = FALSE])
    if (sum(complete_cases) < nrow(data) * 0.5)
      issues$warnings <- c(issues$warnings, glue("Only {sum(complete_cases)} complete cases out of {nrow(data)}"))
  }
  if (length(issues$errors)) { cat("ERRORS:\n"); cat(paste("  -", issues$errors), sep = "\n") }
  if (length(issues$warnings)) { cat("WARNINGS:\n"); cat(paste("  -", issues$warnings), sep = "\n") }
  if (!length(issues$errors) && !length(issues$warnings)) cat("All data validation checks passed\n")
  issues
}


#Prepare data
cat("=== DATA CLEANING AND PREPARATION ===\n")
clean_data <- data %>%
  mutate(across(where(is.character), as.factor)) %>%
  filter(!if_any(any_of(c("biomass_kg_ha","abundance_ind_250m2","FRic","FEve","FDiv","FDis","RaoQ")), is.na)) %>%
  mutate(Site_ID = paste(Fishing_Ground, Dive_Site, sep = "_")) %>%
  group_by(Fishing_Ground) %>% filter(n() >= 3) %>% ungroup()

cat("Data cleaning completed:\n")
cat("  - Original data:", nrow(data), "rows\n")
cat("  - Cleaned data:", nrow(clean_data), "rows\n")
cat("  - Locations:", paste(unique(clean_data$Fishing_Ground), collapse = ", "), "\n")


#Pairwise PERMANOVA

pairwise.adonis2 <- function(x, data, permutations = 999, method = "euclidean", p_adjust = "BH", ...) {
  locations <- unique(data$Fishing_Ground)
  results <- list(); pairwise_results <- data.frame()
  for (i in 1:(length(locations)-1)) for (j in (i+1):length(locations)) {
    loc1 <- locations[i]; loc2 <- locations[j]
    pair_data <- data[data$Fishing_Ground %in% c(loc1, loc2), , drop = FALSE]
    pair_matrix <- x[data$Fishing_Ground %in% c(loc1, loc2), , drop = FALSE]
    if (nrow(pair_matrix) > 2 && ncol(pair_matrix) > 0) {
      dist_matrix <- vegan::vegdist(pair_matrix, method = method)
      permanova <- vegan::adonis2(dist_matrix ~ Fishing_Ground, data = pair_data, permutations = permutations)
      row <- data.frame(
        Comparison = paste(loc1, "vs", loc2),
        F_value = permanova$F[1], R2 = permanova$R2[1],
        p_value = permanova$`Pr(>F)`[1], n = nrow(pair_data)
      )
      pairwise_results <- rbind(pairwise_results, row)
      results[[paste(loc1, "vs", loc2)]] <- list(F_value=row$F_value,R2=row$R2,p_value=row$p_value,n=row$n)
    }
  }
  if (nrow(pairwise_results) > 0) {
    pairwise_results$p_adj <- p.adjust(pairwise_results$p_value, method = p_adjust)
    for (i in 1:nrow(pairwise_results)) results[[pairwise_results$Comparison[i]]]$p_adj <- pairwise_results$p_adj[i]
  }
  list(pairwise_results = pairwise_results, detailed_results = results)
}

#Diagnostics plots
create_model_diagnostics <- function(model, response_var, model_type = "lm") {
  plots <- list(); diagnostics <- list()
  if (inherits(model, "lm") || inherits(model, "lmerMod")) {
    fitted_vals <- fitted(model)
    residuals_vals <- if (inherits(model,"lmerMod")) residuals(model, type="pearson") else residuals(model)
    
    p1 <- ggplot(data.frame(fitted=fitted_vals, residuals=residuals_vals),
                 aes(x=fitted,y=residuals)) +
      geom_point(alpha=.6, size=2,na.rm=TRUE) + geom_hline(yintercept=0, linetype="dashed", color="red") +
      geom_smooth(method="loess", se=TRUE, color="blue", alpha=.2,na.rm=TRUE) +
      labs(title=paste("Residuals vs Fitted -", response_var),
           subtitle=paste("Model type:", model_type), x="Fitted values", y="Residuals") +
      theme_minimal() + theme(plot.title = element_text(face="bold"))
    
    p2 <- ggplot(data.frame(residuals=residuals_vals), aes(sample=residuals)) +
      stat_qq(alpha=.6, size=2) + stat_qq_line(color="red") +
      labs(title=paste("Q-Q Plot -", response_var), subtitle="Checking normality",
           x="Theoretical Quantiles", y="Sample Quantiles") +
      theme_minimal() + theme(plot.title = element_text(face="bold"))
    
    p3 <- ggplot(data.frame(fitted=fitted_vals, sqrt_abs_resid=sqrt(abs(residuals_vals))),
                 aes(x=fitted, y=sqrt_abs_resid)) +
      geom_point(alpha=.6, size=2,na.rm=TRUE) +
      geom_smooth(method="loess", se=TRUE, color="blue", alpha=.2,na.rm=TRUE) +
      labs(title=paste("Scale-Fishing_Ground Plot -", response_var),
           subtitle="Homoscedasticity check", x="Fitted values", y="√|Standardized residuals|") +
      theme_minimal() + theme(plot.title = element_text(face="bold"))
    
    if (inherits(model,"lm")) {
      leverage <- hatvalues(model)
      p4 <- ggplot(data.frame(leverage=leverage, residuals=residuals_vals),
                   aes(x=leverage, y=residuals)) +
        geom_point(alpha=.6, size=2,na.rm=TRUE) + geom_hline(yintercept=0, linetype="dashed", color="red") +
        labs(title=paste("Residuals vs Leverage -", response_var), x="Leverage", y="Residuals") +
        theme_minimal() + theme(plot.title = element_text(face="bold"))
    } else {
      p4 <- ggplot() + theme_void() + labs(title="Leverage plot not available for mixed models")
    }
    plots <- list(residuals_fitted=p1, qq_plot=p2, scale_location=p3, leverage=p4)
    
    if (inherits(model,"lm")) {
      normality_test <- shapiro.test(residuals_vals)
      hetero_test <- tryCatch(car::ncvTest(model), error = function(e) NULL)
      diagnostics$normality_test <- normality_test; diagnostics$heteroscedasticity_test <- hetero_test
      cat("Model diagnostics for", response_var, ":\n")
      cat("  - Shapiro-Wilk W =", round(normality_test$statistic,4),
          "p =", round(normality_test$p.value,4), "\n")
      if (!is.null(hetero_test)) {
        cat("  - NCV Chi² =", round(hetero_test$ChiSquare,4),
            "p =", round(hetero_test$p,4), "\n")
      }
    }
  }
  list(plots=plots, diagnostics=diagnostics)
}

#Partial Plots (Conditional Effects)
#Load required packages for effect plots
if (!require("effects")) {
  install.packages("effects", dependencies = TRUE)
  library(effects)
}
if (!require("ggeffects")) {
  install.packages("ggeffects", dependencies = TRUE)
  library(ggeffects)
}

create_conditional_partial_plots <- function(results, output_prefix = "partial_plots") {
  cat("\n=== CREATING CONDITIONAL PARTIAL PLOTS FOR TOP MODELS ===\n")
  
  #Extract model comparison results
  model_comp <- results$model_comparison
  if (is.null(model_comp) || nrow(model_comp) == 0) {
    cat("No model comparison results found.\n")
    return(NULL)
  }
  
  plot_list <- list()
  
  for (i in 1:nrow(model_comp)) {
    resp <- model_comp$Response[i]
    preferred_model <- model_comp$Preferred[i]
    
    cat("Creating partial plots for:", resp, "- Model:", preferred_model, "\n")
    
    #Get the appropriate model object
    if (preferred_model == "LMM" && !is.null(results$mixed_models[[resp]])) {
      model <- results$mixed_models[[resp]]$model
      model_type <- "LMM"
    } else if (preferred_model == "LM" && !is.null(results$linear_models[[resp]])) {
      model <- results$linear_models[[resp]]$model
      model_type <- "LM"
    } else {
      cat("  Model not found for", resp, "\n")
      next
    }
    
    #Get predictor names from the model
    if (inherits(model, "lmerMod")) {
      predictors <- names(fixef(model))
      predictors <- predictors[!predictors %in% "(Intercept)"]
    } else if (inherits(model, "lm")) {
      predictors <- names(coef(model))
      predictors <- predictors[!predictors %in% "(Intercept)"]
    } else {
      cat("  Unsupported model type for", resp, "\n")
      next
    }
    
    #Clean up predictor names (remove factor level indicators)
    base_predictors <- unique(gsub("Geomorphology[0-9]+", "Geomorphology", 
                                   gsub("crypto5[A-Za-z]+", "crypto", predictors)))
    base_predictors <- base_predictors[!base_predictors %in% c("(Intercept)")]
    
    #Create partial plots for each predictor
    resp_plots <- list()
    
    for (pred in base_predictors) {
      cat("  Plotting:", pred, "\n")
      
      tryCatch({
        if (model_type == "LMM") {
          #For mixed models: conditional effects (including random effects)
          eff_data <- Effect(pred, model)
          eff_df <- as.data.frame(eff_data)
          p_gg <- ggplot(eff_df, aes_string(x = pred, y = "fit")) +
            geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
            geom_line(size = 1) +
            labs(title = paste("Conditional Partial Effect:", pred, "on", resp),
                 subtitle = paste("Model:", model_type),
                 x = pred, y = paste("Predicted", resp)) +
            theme_minimal() +
            theme(plot.title = element_text(face = "bold"))
          
        } else {
          #For linear models: partial residuals plot
          pred_eff <- ggpredict(model, terms = pred)
          p_gg <- plot(pred_eff) +
            labs(title = paste("Partial Effect:", pred, "on", resp),
                 subtitle = paste("Model:", model_type)) +
            theme_minimal() +
            theme(plot.title = element_text(face = "bold"))
        }
        
        #Handle categorical predictors specially
        if (pred %in% c("Geomorphology")) {
          #For categorical variables, use points with error bars
          if (model_type == "LMM") {
            eff_data <- Effect(pred, model)
            eff_df <- as.data.frame(eff_data)
            p_gg <- ggplot(eff_df, aes_string(x = pred, y = "fit")) +
              geom_point(size = 3, color = "blue") +
              geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
              labs(title = paste("Conditional Effect:", pred, "on", resp),
                   subtitle = paste("Model:", model_type, "- Points show estimated marginal means ± CI"),
                   x = pred, y = paste("Predicted", resp)) +
              theme_minimal() +
              theme(plot.title = element_text(face = "bold"),
                    axis.text.x = element_text(angle = 45, hjust = 1))
          }
        }
        
        resp_plots[[pred]] <- p_gg
        print(p_gg)
        
        #Save individual plots
        ggsave(file.path("figures", 
                         paste0(output_prefix, "_", resp, "_", pred, ".png")),
               p_gg, width = 8, height = 6, dpi = 300, bg = "white")
        
      }, error = function(e) {
        cat("    Error creating plot for", pred, ":", conditionMessage(e), "\n")
      })
    }
    
    #Create a combined plot for all predictors if there are multiple
    if (length(resp_plots) > 1) {
      tryCatch({
        combined_plot <- wrap_plots(resp_plots, ncol = 2) + 
          plot_annotation(title = paste("Partial Effects for", resp),
                          subtitle = paste("Top model:", preferred_model))
        
        print(combined_plot)
        ggsave(file.path("figures", 
                         paste0(output_prefix, "_", resp, "_combined.png")),
               combined_plot, width = 12, height = 3 * ceiling(length(resp_plots)/2), 
               dpi = 300, bg = "white")
        
        plot_list[[resp]] <- combined_plot
      }, error = function(e) {
        cat("    Error creating combined plot for", resp, ":", conditionMessage(e), "\n")
      })
    } else if (length(resp_plots) == 1) {
      plot_list[[resp]] <- resp_plots[[1]]
    }
  }
  
  cat("Conditional partial plots completed!\n")
  return(plot_list)
}

# Correlation plot for responses
create_correlation_plot <- function(data, response_vars, method = "pearson") {
  numeric_data <- data %>% dplyr::select(where(is.numeric)) %>% dplyr::select(any_of(response_vars)) %>% na.omit()
  if (ncol(numeric_data) < 2) { cat("Insufficient numeric variables for correlation matrix\n"); return(NULL) }
  cor_matrix <- cor(numeric_data, method = method)
  p <- ggcorrplot::ggcorrplot(cor_matrix, method="circle", type="lower", lab=TRUE, lab_size=3,
                              colors=c("#6D9EC1","white","#E46726"), outline.color="white",
                              show.legend=TRUE, title="Correlation Matrix of Response Variables",
                              ggtheme=theme_minimal()) +
    theme(plot.title = element_text(face="bold", hjust=.5),
          axis.text.x = element_text(angle=45, hjust=1))
  p
}


#Non-parametric tests

nonparametric_analysis <- function(data, response_vars, group_var = "Fishing_Ground") {
  results <- list()
  cat("\n=== NON-PARAMETRIC ANALYSIS ===\n")
  for (resp in response_vars) {
    cat("\n---", resp, "---\n")
    kruskal_formula <- as.formula(paste(resp, "~", group_var))
    kruskal_test <- kruskal.test(kruskal_formula, data = data)
    cat("Kruskal-Wallis: H =", round(kruskal_test$statistic,3),
        "df =", kruskal_test$parameter, "p =", round(kruskal_test$p.value,4), "\n")
    if (kruskal_test$p.value < 0.05) {
      pairwise_wilcox <- pairwise.wilcox.test(data[[resp]], data[[group_var]],
                                              p.adjust.method="BH", exact=FALSE, na.action=na.omit)
      cat("Pairwise Wilcoxon BH-adjusted p-values:\n"); print(pairwise_wilcox$p.value)
    } else pairwise_wilcox <- NULL
    n <- nrow(data); k <- length(unique(data[[group_var]]))
    epsilon_squared <- as.numeric((kruskal_test$statistic - (k - 1)) / (n - k))
    results[[resp]] <- list(kruskal_wallis=kruskal_test, pairwise_wilcoxon=pairwise_wilcox,
                            effect_size=epsilon_squared, n_groups=k, total_n=n)
  }
  results
}


#PCA helper for plots
create_pca_plots_comprehensive <- function(pca_result, metadata, response_vars) {
  pca_scores <- as.data.frame(pca_result$x)
  pca_scores$Fishing_Ground <- if ("Fishing_Ground" %in% colnames(metadata)) metadata$Fishing_Ground else factor("Unknown")
  variance_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
  group_sizes <- table(pca_scores$Fishing_Ground)
  has_ellipse <- names(group_sizes[group_sizes >= 3])
  
  pca_main <- ggplot(pca_scores, aes(x=PC1, y=PC2, color=Fishing_Ground)) +
    geom_point(size=3, alpha=.8) +
    { if (length(has_ellipse) > 0)
      stat_ellipse(data=subset(pca_scores, Fishing_Ground %in% has_ellipse),
                   level=.95, alpha=.2, size=1) else NULL } +
    labs(x=paste0("PC1 (", variance_explained[1], "%)"),
         y=paste0("PC2 (", variance_explained[2], "%)"),
         title="Principal Component Analysis",
         subtitle="Ellipses shown only for groups with n ≥ 3") +
    theme_minimal() + theme(plot.title = element_text(face="bold", hjust=.5)) +
    scale_color_viridis_d(option="plasma")
  
  variance_data <- data.frame(PC=paste0("PC",1:length(variance_explained)),
                              Variance=variance_explained,
                              Cumulative=cumsum(variance_explained)) %>%
    dplyr::slice(1:min(8, n()))
  pca_variance <- ggplot(variance_data, aes(x=reorder(PC,-Variance), y=Variance)) +
    geom_bar(stat="identity", alpha=.7) +
    geom_line(aes(y=Cumulative, group=1), color="red", size=1) +
    geom_point(aes(y=Cumulative), color="red", size=2) +
    labs(title="PCA Variance Explained", subtitle="Bar: Individual PC, Line: Cumulative",
         x="Principal Component", y="Variance Explained (%)") +
    theme_minimal() + theme(plot.title = element_text(face="bold", hjust=.5)) +
    scale_y_continuous(sec.axis = sec_axis(~./100, name="Cumulative Proportion"))
  
  loadings <- as.data.frame(pca_result$rotation[, 1:2, drop=FALSE]); loadings$Variable <- rownames(loadings)
  pca_loadings <- ggplot(loadings, aes(x=PC1, y=PC2)) +
    geom_segment(aes(xend=0, yend=0), arrow=arrow(length=unit(0.2,"cm")), alpha=.7) +
    geom_text(aes(label=Variable), hjust=-.1, vjust=-.1, size=3, fontface="bold") +
    labs(title="PCA Variable Loadings",
         x=paste0("PC1 (", variance_explained[1], "%)"),
         y=paste0("PC2 (", variance_explained[2], "%)")) +
    theme_minimal() + theme(plot.title = element_text(face="bold", hjust=.5))
  
  list(main=pca_main, variance=pca_variance, loadings=pca_loadings)
}

#Collinearity control + scaling for LM & LMM
#1)Scale numeric predictors in-place (responses are NOT scaled)
scale_predictors_inplace <- function(df, predictor_vars, center = TRUE, scale = TRUE) {
  num_preds <- predictor_vars[predictor_vars %in% names(df)]
  num_preds <- num_preds[sapply(num_preds, \(v) is.numeric(df[[v]]))]
  if (!length(num_preds)) return(list(data=df, scaled_cols=character(0)))
  df[num_preds] <- lapply(df[num_preds], \(z) as.numeric(scale(z, center=center, scale=scale)))
  list(data=df, scaled_cols=num_preds)
}

#2) Correlation screen (numeric-only): remove variables to keep r <= threshold
correlation_screen <- function(df, predictors, thr = 0.75) {
  keep <- predictors
  numeric_preds <- keep[keep %in% names(df)]
  numeric_preds <- numeric_preds[sapply(numeric_preds, \(v) is.numeric(df[[v]]))]
  if (length(numeric_preds) < 2) return(list(kept=keep, dropped_pairs=data.frame()))
  cm <- stats::cor(df[, numeric_preds, drop=FALSE], use="pairwise.complete.obs")
  dropped <- data.frame(var_removed=character(), rival=character(), r=numeric())
  #Greedy: drop the variable with the largest average r when any pair exceeds threshold
  repeat {
    above <- which(abs(cm) > thr & lower.tri(cm), arr.ind = TRUE)
    if (!nrow(above)) break
    pairs <- tibble(v1=rownames(cm)[above[,1]], v2=colnames(cm)[above[,2]],
                    r = mapply(\(i,j) cm[i,j], above[,1], above[,2]))
    involved <- unique(c(pairs$v1, pairs$v2))
    m <- abs(cm[involved, involved, drop=FALSE]); diag(m) <- NA_real_
    avg_abs <- sort(colMeans(m, na.rm = TRUE), decreasing = TRUE)
    drop_var <- names(avg_abs)[1]
    rival_row <- pairs %>% filter(v1 == drop_var | v2 == drop_var) %>% slice_max(abs(r), n=1)
    dropped <- rbind(dropped, data.frame(var_removed=drop_var,
                                         rival=ifelse(rival_row$v1==drop_var, rival_row$v2, rival_row$v1),
                                         r=rival_row$r))
    keep <- setdiff(keep, drop_var)
    numeric_preds <- setdiff(numeric_preds, drop_var)
    if (length(numeric_preds) < 2) break
    cm <- stats::cor(df[, numeric_preds, drop=FALSE], use="pairwise.complete.obs")
  }
  list(kept=keep, dropped_pairs=dropped)
}

#3)Iterative (G)VIF screen using LM (proxy for LMM fixedpart collinearity)
vif_screen <- function(df, response, predictors, vif_threshold = 5, exclude_terms = character(0)) {
  preds <- setdiff(predictors, exclude_terms)
  dropped <- tibble::tibble(term=character(), GVIF=numeric(), Df=integer(), GVIF_adj=numeric())
  path <- list()
  
  repeat {
    if (!length(preds)) break
    form <- as.formula(paste(response, "~", paste(preds, collapse=" + ")))
    fit  <- tryCatch(lm(form, data=df), error=function(e) NULL)
    if (is.null(fit)) break
    
    vf <- tryCatch(car::vif(fit), error=function(e) NULL)
    if (is.null(vf)) break
    
    vif_tbl <- if (is.matrix(vf)) {
      out <- as.data.frame(vf)
      out$Term <- rownames(out)
      out$GVIF_adj <- out$GVIF^(1/(2*out$Df))
      rownames(out) <- NULL
      out[, c("Term","Df","GVIF","GVIF_adj")]
    } else {
      data.frame(Term=names(vf), Df=1L, GVIF=as.numeric(vf), GVIF_adj=sqrt(as.numeric(vf)))
    }
    
    #thresholding rule applied here
    vif_tbl$Check <- ifelse(vif_tbl$Df == 1, vif_tbl$GVIF, vif_tbl$GVIF_adj)
    path[[length(path)+1]] <- vif_tbl
    
    worst <- vif_tbl[which.max(vif_tbl$Check), ]
    if (worst$Check <= vif_threshold) break
    
    preds <- setdiff(preds, worst$Term)
    dropped <- dplyr::bind_rows(dropped, worst[, c("Term","GVIF","Df","GVIF_adj")])
  }
  
  list(kept=preds, dropped_table=dropped, path=path)
}


#p.value for lmer models
get_lmm_tidy_with_p <- function(mod, effect = "fixed") {
  td <- tryCatch(broom.mixed::tidy(mod, effects = effect), error = function(e) NULL)
  if (is.null(td)) return(NULL)
  if (!"p.value" %in% names(td) || all(is.na(td$p.value))) {
    cs <- tryCatch(coef(summary(mod)), error = function(e) NULL)
    if (!is.null(cs)) {
      pcol <- intersect(colnames(cs), c("Pr(>|t|)", "p.value", "pval", "p"))
      if (length(pcol) == 1) {
        p_map <- tibble::tibble(term = rownames(cs), p.value = as.numeric(cs[, pcol]))
        td <- dplyr::left_join(td, p_map, by = "term")
      } else {
        td$p.value <- NA_real_
      }
    } else {
      td$p.value <- NA_real_
    }
  }
  td
}

compute_lm_robust <- function(lm_model) {
  vc <- sandwich::vcovHC(lm_model, type = "HC3")
  ct <- lmtest::coeftest(lm_model, vcov. = vc)
  tb <- broom::tidy(ct)
  names(tb)[names(tb)=="statistic"] <- "t.value"
  names(tb)[names(tb)=="std.error"] <- "robust.se"
  tb
}


#Univariate panel (LM + LMM) for excluded (or all) predictors

run_univariate_panel <- function(analysis_data,
                                 response_vars,
                                 all_predictors,
                                 final_fixed_predictors,
                                 random_effects = c("(1|Dive_Site)"),
                                 output_prefix = "univariate",
                                 scope = c("excluded_only","all"),
                                 save_results = TRUE) {
  scope <- match.arg(scope)
  preds <- if (scope == "excluded_only") setdiff(all_predictors, final_fixed_predictors) else all_predictors
  preds <- unique(preds[preds %in% names(analysis_data)])
  preds <- preds[!preds %in% c("Fishing_Ground")]
  if (!length(preds)) { message("Univariate panel: nothing to test."); return(invisible(NULL)) }
  
  out_lm <- list(); out_lmm <- list()
  
  for (resp in response_vars) {
    if (!resp %in% names(analysis_data)) next
    for (pr in preds) {
      df <- analysis_data[, c(resp, pr, "Fishing_Ground"), drop = FALSE]
      df <- df[stats::complete.cases(df), , drop = FALSE]
      if (!nrow(df)) next
      
      #LM
      f_lm <- stats::as.formula(paste(resp, "~", pr))
      fit_lm <- tryCatch(stats::lm(f_lm, data = df), error = function(e) NULL)
      if (!is.null(fit_lm)) {
        sm <- broom::tidy(fit_lm)
        gl <- broom::glance(fit_lm)
        sm$Response <- resp; sm$Predictor <- pr; gl$Response <- resp; gl$Predictor <- pr
        out_lm[[paste(resp, pr, "lm", sep = "_")]] <- list(tidy = sm, glance = gl)
      }
      
      #LMM with multiple locations 
      if ("Fishing_Ground" %in% names(df) && length(unique(df$Fishing_Ground)) > 1) {
        f_lmm <- stats::as.formula(paste(resp, "~", pr, "+ (1|Dive_Site)"))
        fit_lmm <- tryCatch(
          lmerTest::lmer(f_lmm, data = df,
                         control = lme4::lmerControl(optimizer = "bobyqa",
                                                     optCtrl = list(maxfun = 100000))),
          error = function(e) NULL
        )
        if (!is.null(fit_lmm)) {
          sm <- get_lmm_tidy_with_p(fit_lmm, effect = "fixed")
          gl <- tryCatch(broom.mixed::glance(fit_lmm), error = function(e) NULL)
          if (!is.null(sm)) {
            sm$Response <- resp; sm$Predictor <- pr
            out_lmm[[paste(resp, pr, "lmm", sep = "_")]] <- list(tidy = sm, glance = gl)
          }
        }
      }
    }
  }
  
  bind_and_adjust <- function(lst) {
    if (!length(lst)) return(NULL)
    td <- purrr::map_dfr(lst, "tidy")
    gl <- purrr::map_dfr(lst, "glance")
    td_pred <- dplyr::semi_join(td, dplyr::distinct(td, Response, Predictor),
                                by = c("Response", "term" = "Predictor"))
    if (!"p.value" %in% names(td_pred)) td_pred$p.value <- NA_real_
    if (nrow(td_pred)) {
      td_pred <- td_pred %>%
        dplyr::group_by(Response) %>%
        dplyr::mutate(p_adj_BH = if (all(is.na(p.value))) NA_real_
                      else stats::p.adjust(replace(p.value, is.na(p.value), 1), method = "BH"),
                      sig = dplyr::case_when(
                        !is.na(p_adj_BH) & p_adj_BH < 0.001 ~ "***",
                        !is.na(p_adj_BH) & p_adj_BH < 0.01  ~ "**",
                        !is.na(p_adj_BH) & p_adj_BH < 0.05  ~ "*",
                        !is.na(p_adj_BH) & p_adj_BH < 0.1   ~ ".",
                        TRUE ~ ""
                      )) %>%
        dplyr::ungroup()
    }
    list(tidy = td, tidy_predictor = td_pred, glance = gl)
  }
  
  lm_res  <- bind_and_adjust(out_lm)
  lmm_res <- bind_and_adjust(out_lmm)
  
  if (isTRUE(save_results)) {
    if (!is.null(lm_res)) {
      if (!is.null(lm_res$tidy) && nrow(lm_res$tidy))
        utils::write.csv(lm_res$tidy, file.path("tables", paste0(output_prefix, "_LM_all_terms.csv")), row.names = FALSE)
      if (!is.null(lm_res$tidy_predictor) && nrow(lm_res$tidy_predictor))
        utils::write.csv(lm_res$tidy_predictor, file.path("tables", paste0(output_prefix, "_LM_predictor_only_BH.csv")), row.names = FALSE)
      if (!is.null(lm_res$glance) && nrow(lm_res$glance))
        utils::write.csv(lm_res$glance, file.path("tables", paste0(output_prefix, "_LM_glance.csv")), row.names = FALSE)
    }
    if (!is.null(lmm_res)) {
      if (!is.null(lmm_res$tidy) && nrow(lmm_res$tidy))
        utils::write.csv(lmm_res$tidy, file.path("tables", paste0(output_prefix, "_LMM_all_terms.csv")), row.names = FALSE)
      if (!is.null(lmm_res$tidy_predictor) && nrow(lmm_res$tidy_predictor))
        utils::write.csv(lmm_res$tidy_predictor, file.path("tables", paste0(output_prefix, "_LMM_predictor_only_BH.csv")), row.names = FALSE)
      if (!is.null(lmm_res$glance) && nrow(lmm_res$glance))
        utils::write.csv(lmm_res$glance, file.path("tables", paste0(output_prefix, "_LMM_glance.csv")), row.names = FALSE)
    }
  }
  
  invisible(list(lm = lm_res, lmm = lmm_res))
}


#CONFIGURATIONSu 
##this is where you modify your predictir vars
analysis_config <- list(
  core_ecological = list(
    response_vars = c("biomass_kg_ha", "abundance_ind_250m2", "FRic", "FEve"),
    predictor_vars = c("Geomorphology","sedimnt","nutrint","NO_TK_AREA","Fishing_Ground",
                       "crypto5BROF","crypto5BRIF","crypto5LR","crypto5BRin","crypto5BRout",
                       "sst_sd_6_year_mean","sst_q90_6_year_mean","sst_mean_6_year_mean"),
    random_effects = c("(1|Dive_Site)"),
    description = "Core ecological variables analysis"
  ),
  functional_diversity = list(
    response_vars = c("FRic","FEve","FDiv","FDis","RaoQ"),
    predictor_vars = c("Geomorphology","sedimnt","nutrint","NO_TK_AREA","Fishing_Ground",
                       "crypto5BROF","crypto5BRIF","crypto5LR","crypto5BRin","crypto5BRout",
                       "sst_sd_6_year_mean","sst_q90_6_year_mean","sst_mean_6_year_mean"),
    random_effects = c("(1|Dive_Site)"),
    description = "Functional diversity metrics analysis"
  ),
  biomass_focused = list(
    response_vars = c("biomass_kg_ha","abundance_ind_250m2"),
    predictor_vars = c("Geomorphology","sedimnt","nutrint","NO_TK_AREA","Fishing_Ground",
                       "crypto5BROF","crypto5BRIF","crypto5LR","crypto5BRin","crypto5BRout",
                       "sst_sd_6_year_mean","sst_q90_6_year_mean","sst_mean_6_year_mean"),
    random_effects = c("(1|Dive_Site)"),
    description = "Biomass and abundance focused analysis"
  ),
  all_responses = list(
    response_vars = c("biomass_kg_ha","abundance_ind_250m2","FRic","FEve","FDiv","FDis","RaoQ"),
    predictor_vars = c("Geomorphology","sedimnt","nutrint","NO_TK_AREA","Fishing_Ground",
                       "crypto5BROF","crypto5BRIF","crypto5LR","crypto5BRin","crypto5BRout",
                       "sst_sd_6_year_mean","sst_q90_6_year_mean","sst_mean_6_year_mean"),
    random_effects = c("(1|Dive_Site)"),
    description = "all response variables"
  )
)

#MAIN ANALYSIS FUNCTION
analyze_reef_data_enhanced <- function(
    response_vars,
    predictor_vars,
    random_effects = c("(1|Dive_Site)"),
    data = clean_data,
    output_prefix = "_analysis",
    create_plots = TRUE,
    save_results = TRUE,
    run_diagnostics = TRUE,
    run_nonparametric = TRUE,
    # model-build options
    lmm_location_as_random_only = TRUE,   # drop 'Fishing_Ground' from fixed if used as random
    scale_predictors_for_LMM = TRUE,      # ALWAYS recommended
    scale_predictors_for_LM  = TRUE,      # generally recommended
    cor_threshold = 0.75,                 # numeric pairwise r limit
    vif_threshold = 5,                    # GVIF^(1/(2*Df)) limit (≈2)
    run_univariate_scope = "excluded_only", # "excluded_only" or "all"
    create_conditional_partial_plots = TRUE
) {
  cat("\n", strrep("=", 74), "\n")
  cat("STARTING ANALYSIS\n")
  cat(strrep("=", 74), "\n")
  
  if (length(response_vars) == 0) stop("No response variables specified.")
  if (length(predictor_vars) == 0) stop("No predictor variables specified.")
  
  #Validation 
  validation_issues <- validate_analysis_data(data, response_vars, predictor_vars)
  if (length(validation_issues$errors) > 0) stop("Critical validation errors. Aborting.")
  
  #Build analysis dataset
  analysis_data <- data %>%
    dplyr::select(dplyr::any_of(c(response_vars, predictor_vars, "Site_ID","Transect","Dive_Site","Fishing_Ground"))) %>%
    dplyr::filter(stats::complete.cases(.))
  if (nrow(analysis_data) == 0) stop("No complete cases after filtering. Check inputs.")
  
  #Scaling (predictors only). LMM always; LM optional; we scale once here for both.
  to_scale <- predictor_vars[predictor_vars %in% names(analysis_data)]
  to_scale <- to_scale[sapply(to_scale, \(v) is.numeric(analysis_data[[v]]))]
  do_scale <- (scale_predictors_for_LMM || scale_predictors_for_LM)
  scaled_cols <- character(0)
  if (do_scale && length(to_scale)) {
    sc <- scale_predictors_inplace(analysis_data, predictor_vars = to_scale, center = TRUE, scale = TRUE)
    analysis_data <- sc$data; scaled_cols <- sc$scaled_cols
  }
  
  cat("Analysis dataset created:\n")
  cat("  - Variables:", length(unique(c(response_vars, predictor_vars))), "\n")
  cat("  - Observations:", nrow(analysis_data), "\n")
  cat("  - Locations:", paste(utils::head(unique(analysis_data$Fishing_Ground), 20), collapse = ", "),
      if (length(unique(analysis_data$Fishing_Ground)) > 20) " ... [truncated]" else "", "\n", sep = "")
  if (length(scaled_cols)) cat("  - Scaled predictors:", paste(scaled_cols, collapse = ", "), "\n")
  
  #Store meta data
  results <- list(metadata = list(
    response_vars = response_vars, predictor_vars = predictor_vars, random_effects = random_effects,
    sample_size = nrow(analysis_data), locations = unique(analysis_data$Fishing_Ground),
    n_locations = length(unique(analysis_data$Fishing_Ground)), timestamp = Sys.time(),
    validation_issues = validation_issues, scaled_predictors = scaled_cols
  ))
  
  #Descriptive  stats
  cat("\n=== DESCRIPTIVE STATISTICS ===\n")
  location_summary <- analysis_data %>%
    group_by(.data$Fishing_Ground) %>%
    summarise(
      n_transects = n(),
      across(all_of(response_vars),
             list(mean=~mean(.x, na.rm=TRUE),
                  sd=~sd(.x, na.rm=TRUE),
                  median=~median(.x, na.rm=TRUE),
                  min=~min(.x, na.rm=TRUE),
                  max=~max(.x, na.rm=TRUE),
                  cv=~ifelse(mean(.x, na.rm=TRUE)==0, NA_real_, sd(.x, na.rm=TRUE)/mean(.x, na.rm=TRUE))),
             .names = "{.col}_{.fn}"),
      .groups="drop"
    )
  overall_summary <- analysis_data %>%
    summarise(
      n_total = n(), n_locations = n_distinct(.data$Fishing_Ground),
      across(all_of(response_vars),
             list(mean=~mean(.x, na.rm=TRUE),
                  sd=~sd(.x, na.rm=TRUE),
                  median=~median(.x, na.rm=TRUE),
                  min=~min(.x, na.rm=TRUE),
                  max=~max(.x, na.rm=TRUE),
                  cv=~ifelse(mean(.x, na.rm=TRUE)==0, NA_real_, sd(.x, na.rm=TRUE)/mean(.x, na.rm=TRUE))),
             .names = "{.col}_{.fn}")
    )
  print(location_summary)
  results$descriptive <- list(by_location = location_summary, overall = overall_summary)
  if (save_results) {
    write.csv(location_summary, file.path("tables", paste0(output_prefix, "_location_summary.csv")), row.names = FALSE)
    write.csv(overall_summary,   file.path("tables", paste0(output_prefix, "_overall_summary.csv")),   row.names = FALSE)
  }
  
  #Visualisation 
  plots <- list()
  if (create_plots) {
    cat("\n=== VISUALIZATION ===\n")
    for (resp in response_vars) if (resp %in% names(analysis_data)) {
      p_box <- ggplot(analysis_data, aes(x=Fishing_Ground, y=.data[[resp]], fill=Fishing_Ground)) +
        geom_boxplot(alpha=.7, outlier.shape=NA) +
        geom_jitter(width=.2, alpha=.6, size=1.5) +
        stat_summary(fun=mean, geom="point", shape=18, size=3, color="red") +
        labs(title=paste("Distribution of", resp, "by Fishing_Ground"),
             subtitle=glue("n = {nrow(analysis_data)} transects; red diamond = mean"),
             y=resp, x="Fishing_Ground") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
        scale_fill_viridis_d(option="viridis")
      plots[[paste0("boxplot_",resp)]] <- p_box; print(p_box)
      if (save_results) ggsave(file.path("figures", paste0(output_prefix, "_boxplot_", resp, ".png")),
                               p_box, width=10, height=6, dpi=300, bg="white")
    }
    cor_plot <- create_correlation_plot(analysis_data, response_vars)
    if (!is.null(cor_plot)) {
      plots$correlation_matrix <- cor_plot; print(cor_plot)
      if (save_results) ggsave(file.path("figures", paste0(output_prefix, "_correlation_matrix.png")),
                               cor_plot, width=8, height=6, dpi=300, bg="white")
    }
  }
  results$plots <- plots
  
  #Non-parametric 
  if (run_nonparametric) {
    nonparametric_results <- nonparametric_analysis(analysis_data, response_vars)
    results$nonparametric <- nonparametric_results
    if (save_results) {
      nonparametric_df <- purrr::map_dfr(names(nonparametric_results), function(resp) {
        res <- nonparametric_results[[resp]]
        data.frame(Response=resp, Test="Kruskal-Wallis",
                   Statistic=as.numeric(res$kruskal_wallis$statistic),
                   df=as.numeric(res$kruskal_wallis$parameter),
                   p_value=as.numeric(res$kruskal_wallis$p.value),
                   Effect_Size=as.numeric(res$effect_size),
                   n_groups=res$n_groups, total_n=res$total_n)
      })
      write.csv(nonparametric_df, file.path("tables", paste0(output_prefix, "_nonparametric_results.csv")), row.names = FALSE)
    }
  }
  

  #Collinearity-aware fixed-effect set for LM & LMM
  fixed_candidates <- predictor_vars
  if (lmm_location_as_random_only) fixed_candidates <- setdiff(fixed_candidates, "Fishing_Ground")
  
  #A:Correlation screen
  cor_scr <- correlation_screen(analysis_data, predictors = fixed_candidates, thr = cor_threshold)
  fixed_after_cor <- cor_scr$kept
  
  #B:VIF screen (use first available numeric response as proxy)
  proxy_resp <- response_vars[response_vars %in% names(analysis_data)][1]
  vif_scr <- vif_screen(analysis_data, response = proxy_resp, predictors = fixed_after_cor,
                        vif_threshold = vif_threshold, exclude_terms = character(0))
  fixed_screened <- vif_scr$kept
  
  #Log screening
  screening_log <- list(
    fixed_initial = fixed_candidates,
    correlation_dropped = cor_scr$dropped_pairs,
    vif_dropped = vif_scr$dropped_table,
    final_fixed = fixed_screened
  )
  if (save_results) {
    if (nrow(cor_scr$dropped_pairs)) {
      write.csv(cor_scr$dropped_pairs,
                file.path("tables", paste0(output_prefix, "_COR_dropped_pairs.csv")), row.names = FALSE)
    }
    if (nrow(vif_scr$dropped_table)) {
      write.csv(vif_scr$dropped_table,
                file.path("tables", paste0(output_prefix, "_VIF_dropped_terms.csv")), row.names = FALSE)
    }
    if (length(vif_scr$path)) {
      for (i in seq_along(vif_scr$path)) {
        vt <- vif_scr$path[[i]]
        if (!is.null(vt)) {
          write.csv(vt, file.path("tables", paste0(output_prefix, "_VIF_path_step_", sprintf("%02d", i), ".csv")), row.names = FALSE)
        }
      }
    }
    writeLines(paste("Final fixed predictors:", paste(fixed_screened, collapse = ", ")),
               con = file.path("tables", paste0(output_prefix, "_final_fixed_predictors.txt")))
  }
  
  if (!length(fixed_screened)) {
    cat("⚠️  No fixed-effect predictors remain after collinearity screening. Models will be intercept-only.\n")
  } else {
    cat("Fixed-effect predictors after screening:\n  - ", paste(fixed_screened, collapse = ", "), "\n", sep = "")
  }
  
  #Linear Models (LM)
  cat("\n=== LINEAR MODELS (with diagnostics & VIF) ===\n")
  lm_results <- list(); lm_diagnostics <- list()
  for (resp in response_vars) if (resp %in% names(analysis_data)) {
    cat("\n--- LM for:", resp, "---\n")
    fixed_terms_lm <- fixed_screened  
    form_lm <- as.formula(paste(resp, "~", if (length(fixed_terms_lm)) paste(fixed_terms_lm, collapse=" + ") else "1"))
    lm_model <- lm(form_lm, data = analysis_data)
    lm_summary <- summary(lm_model); print(lm_summary)
    
    anova_result <- tryCatch(car::Anova(lm_model, type = "III"),
                             error=function(e) tryCatch(car::Anova(lm_model, type="II"), error=function(e2) anova(lm_model)))
    if (!is.null(anova_result)) { cat("\nANOVA Results:\n"); print(anova_result) }
    
    if (run_diagnostics) {
      d <- create_model_diagnostics(lm_model, resp, "LM"); lm_diagnostics[[resp]] <- d
      if (save_results && create_plots) {
        for (nm in names(d$plots)) {
          ggsave(file.path("diagnostics", paste0(output_prefix, "_lm_", resp, "_", nm, ".png")),
                 d$plots[[nm]], width=8, height=6, dpi=300, bg="white")
        }
      }
    }
    
    #VIF 
    vif_tbl <- NULL

    vif_tbl <- tryCatch({
      v <- car::vif(lm_model)
      if (is.matrix(v)) {
        out <- data.frame(
          Term     = rownames(v),
          Df       = v[, "Df"],
          GVIF     = v[, "GVIF"],
          GVIF_adj = v[, "GVIF^(1/(2*Df))"],
          row.names = NULL
        )
      } else {
        out <- data.frame(
          Term     = names(v),
          Df       = 1L,
          GVIF     = as.numeric(v),
          GVIF_adj = as.numeric(v),
          row.names = NULL
        )
      }
      out$Check <- ifelse(out$Df > 1, out$GVIF_adj, out$GVIF)
      out$Flag  <- out$Check > 5
      out
    }, error = function(e) {
      message("VIF failed: ", conditionMessage(e))
      data.frame(Term=character(), Df=integer(), GVIF=double(),
                 GVIF_adj=double(), Check=double(), Flag=logical())
    })
    
    cat("\nVIF summary for ", resp, ":\n", sep = "")
    print(vif_tbl)

    
    lm_results[[resp]] <- list(
      model = lm_model,
      summary = lm_summary,
      anova = anova_result,
      coefficients = broom::tidy(lm_model),
      performance  = broom::glance(lm_model),
      diagnostics  = if (run_diagnostics) lm_diagnostics[[resp]] else NULL,
      vif = vif_tbl
    )
    
    robust_tbl <- compute_lm_robust(lm_model)
    if (save_results) {
      utils::write.csv(robust_tbl, file.path("tables", paste0(output_prefix, "_lm_", resp, "_HC3.csv")), row.names = FALSE)
    }
    lm_results[[resp]]$robust_HC3 <- robust_tbl
  }
  results$linear_models <- lm_results
  if (save_results && length(lm_results)) {
    lm_coefs <- purrr::map_dfr(lm_results, \(x) x$coefficients, .id="Response")
    write.csv(lm_coefs, file.path("tables", paste0(output_prefix, "_linear_model_coefficients.csv")), row.names = FALSE)
    lm_perf  <- purrr::map_dfr(lm_results, \(x) x$performance, .id="Response")
    write.csv(lm_perf,  file.path("tables", paste0(output_prefix, "_linear_model_performance.csv")),   row.names = FALSE)
    all_vif <- purrr::map_dfr(names(lm_results), \(nm) { vt <- lm_results[[nm]]$vif; if (is.null(vt)) return(NULL); cbind(Response=nm, vt) })
    if (nrow(all_vif) > 0) write.csv(all_vif, file.path("tables", paste0(output_prefix, "_VIF_all_responses.csv")), row.names = FALSE)
  }
  
  #Mixed Models (LMM) 
  cat("\n=== MIXED MODELS (LMM) WITH DIAGNOSTICS ===\n")
  mixed_results <- list(); mixed_diagnostics <- list()
  singular_flags <- setNames(rep(NA, length(response_vars)), response_vars)
  
  for (resp in response_vars) if (resp %in% names(analysis_data)) {
    cat("\n--- LMM for:", resp, "---\n")
    fixed_terms <- fixed_screened
    fixed_terms <- if (length(fixed_terms)) paste(fixed_terms, collapse=" + ") else "1"
    mixed_formula <- as.formula(paste(resp, "~", fixed_terms, "+", paste(random_effects, collapse=" + ")))
    mixed_model <- tryCatch(
      lmerTest::lmer(
        mixed_formula,
        data = analysis_data,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
      ),
      error = function(e) {
        cat("Complex LMM failed, trying simpler structure with (1|Dive_Site):\n")
        simple_formula <- as.formula(paste(resp, "~", fixed_terms, "+ (1|Dive_Site)"))
        lmerTest::lmer(
          simple_formula,
          data = analysis_data,
          control = lme4::lmerControl(optimizer = "bobyqa")
        )
      }
    )
    
    if (!inherits(mixed_model, "try-error") && !is.null(mixed_model)) {
      mixed_summary <- summary(mixed_model); print(mixed_summary)
      r2 <- tryCatch(MuMIn::r.squaredGLMM(mixed_model), error=function(e) NULL)
      if (!is.null(r2)) { cat("\nR-squared (Marginal / Conditional):\n"); print(r2) }
      varcomp <- as.data.frame(VarCorr(mixed_model)); cat("\nVariance Components:\n"); print(varcomp)
      
      #Type III tests: prefer Kenward–Roger, fall back to Satterthwaite
      primary_ddf <- if (requireNamespace("pbkrtest", quietly = TRUE)) {
        "Kenward-Roger"
      } else {
        "Satterthwaite"
      }
      
      attempt_type3 <- function(ddf) {
        tryCatch(
          anova(mixed_model, type = 3, ddf = ddf),  # bare anova() for S3 dispatch via lmerTest
          error = function(e) { 
            message("Type III (", ddf, ") failed: ", conditionMessage(e))
            NULL
          }
        )
      }
      
      type3_tab   <- attempt_type3(primary_ddf)
      used_method <- primary_ddf
      
      if (is.null(type3_tab) && identical(primary_ddf, "Kenward-Roger")) {
        message("Type III (Kenward–Roger) failed; falling back to Satterthwaite.")
        type3_tab   <- attempt_type3("Satterthwaite")
        used_method <- "Satterthwaite"
      }
      
      if (is.null(type3_tab)) {
        message("Type III failed for both Kenward–Roger and Satterthwaite.")
      } else {
        cat("\nType III tests (", used_method, "):\n", sep = ""); print(type3_tab)
        if (isTRUE(save_results)) {
          utils::write.csv(
            as.data.frame(type3_tab),
            file = file.path(
              "tables",
              sprintf("%s_lmer_%s_TypeIII_%s.csv",
                      output_prefix, resp, gsub("\\s+","", used_method))
            ),
            row.names = TRUE
          )
        }
      }
      
      
      #Auto-pick model & label the test (place right after the KR Type III table)
      AIC_lm  <- suppressWarnings(tryCatch(AIC(lm_model),  error = function(e) NA_real_))
      AIC_lmm <- suppressWarnings(tryCatch(AIC(mixed_model), error = function(e) NA_real_))
      sg      <- suppressWarnings(tryCatch(lme4::isSingular(mixed_model, tol = 1e-5), error = function(e) NA))
      
      delta <- AIC_lmm - AIC_lm
      pick_model <- if (isTRUE(sg) || is.na(AIC_lmm) || (!is.na(delta) && delta >= 2)) {
        "LM"
      } else if (!isTRUE(sg) && !is.na(delta) && delta <= -2) {
        "LMM"
      } else {
        "LM"
      }
      
      cat(sprintf("\nModel choice for %s: %s (ΔAIC = %.2f; singular = %s)\n",
                  resp, pick_model, delta, ifelse(is.na(sg), "NA", as.character(sg))))
      
      choice_tbl <- rbind(choice_tbl,
                          data.frame(response = resp,
                                     AIC_LM = AIC_lm, AIC_LMM = AIC_lmm,
                                     delta = delta,
                                     singular = ifelse(is.na(sg), FALSE, sg),
                                     choice = pick_model,
                                     stringsAsFactors = FALSE)
      )
     
      # Mark singularfits
      singular_flags[resp] <- tryCatch(lme4::isSingular(mixed_model, tol = 1e-5), error = function(e) NA)
      if (isTRUE(singular_flags[resp])) cat("⚠️  Singular LMM for ", resp, ": random-effect variance ~ 0; LM likely preferable.\n", sep="")
      
      if (run_diagnostics && !isTRUE(singular_flags[resp])) {
        dh <- tryCatch(DHARMa::simulateResiduals(mixed_model, n = 1000), error=function(e) NULL)
        if (!is.null(dh)) {
          mixed_diagnostics[[resp]] <- list(dharma = dh)
          if (save_results && create_plots) {
            png(file.path("diagnostics", paste0(output_prefix, "_lmer_", resp, "_DHARMa.png")),
                width=900, height=900, res=150); plot(dh); dev.off()
          }
        }
      }
      
      mixed_results[[resp]] <- list(
        model = mixed_model,
        summary = mixed_summary,
        r2 = r2,
        varcomp = varcomp,
        type3 = type3_tab,
        diagnostics = if (isTRUE(run_diagnostics)) mixed_diagnostics[[resp]] else NULL
      )
      
    } else cat("LMM failed for response:", resp, "\n")
  }
  results$mixed_models <- mixed_results
  if (save_results && length(mixed_results)) {
    sink(file.path("tables", paste0(output_prefix, "_mixed_model_summaries.txt")))
    for (nm in names(mixed_results)) {
      cat("\n====== Mixed Model:", nm, "======\n")
      print(mixed_results[[nm]]$summary)
      if (!is.null(mixed_results[[nm]]$r2)) { cat("\nR2 (marginal/conditional):\n"); print(mixed_results[[nm]]$r2) }
      cat("\nVariance Components:\n"); print(mixed_results[[nm]]$varcomp)
      cat("\nType III tests (Satterthwaite):\n"); print(mixed_results[[nm]]$type3)
    }
    sink()
    lapply(names(mixed_results), function(nm)
      saveRDS(mixed_results[[nm]]$model, file.path("models", paste0(output_prefix, "_lmer_", nm, ".rds"))))
  }
  
  #Model comparison: LM vs LMM (AIC/BIC)
  cat("\n=== MODEL COMPARISON (AIC/BIC): LM vs LMM ===\n")
  model_comp_rows <- list()
  
  if (length(lm_results) > 0 && length(mixed_results) > 0) {
    for (resp in names(lm_results)) {
      lm_obj  <- lm_results[[resp]]$model
      lmm_obj <- if (!is.null(mixed_results[[resp]])) mixed_results[[resp]]$model else NULL
      if (is.null(lm_obj) || is.null(lmm_obj)) { 
        cat("No comparable pair (LM/LMM) for:", resp, "\n"); next 
      }
      
      is_sing <- tryCatch(lme4::isSingular(lmm_obj, tol=1e-5), error=function(e) NA)
      
      # ML refit for the LMM for fair AIC/BIC
      lmm_ml <- tryCatch(lme4::refitML(lmm_obj), error=function(e) {
        tryCatch(update(lmm_obj, REML = FALSE), error=function(e2) NULL)
      })
      
      aic_lm <- AIC(lm_obj); bic_lm <- BIC(lm_obj)
      if (isTRUE(is_sing) || is.null(lmm_ml)) {
        aic_lmm <- NA_real_; bic_lmm <- NA_real_
        delta_aic <- NA_real_; delta_bic <- NA_real_
      } else {
        aic_lmm <- AIC(lmm_ml); bic_lmm <- BIC(lmm_ml)
        delta_aic <- aic_lmm - aic_lm
        delta_bic <- bic_lmm - bic_lm
      }
      
      cat(sprintf("  %s: LM AIC=%.2f, LMM AIC=%s | LM BIC=%.2f, LMM BIC=%s\n",
                  resp, aic_lm,
                  ifelse(is.na(aic_lmm), "NA", sprintf("%.2f", aic_lmm)),
                  bic_lm,
                  ifelse(is.na(bic_lmm), "NA", sprintf("%.2f", bic_lmm))))
      
      model_comp_rows[[resp]] <- data.frame(
        Response  = resp,
        LM_AIC    = aic_lm,  LMM_AIC = aic_lmm,  Delta_AIC = delta_aic,
        LM_BIC    = bic_lm,  LMM_BIC = bic_lmm,  Delta_BIC = delta_bic,
        Singular_LMM = isTRUE(is_sing),
        Preferred = dplyr::case_when(
          isTRUE(is_sing)            ~ "LM",
          is.na(delta_aic)           ~ "LM",
          delta_aic < -2             ~ "LMM",
          delta_aic >  2             ~ "LM",
          TRUE                       ~ "Tie"
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  
  model_comp_df <- dplyr::bind_rows(model_comp_rows)
  results[["model_comparison"]] <- model_comp_df
  if (isTRUE(save_results) && nrow(model_comp_df) > 0) {
    utils::write.csv(model_comp_df, file = file.path("tables", paste0(output_prefix, "_model_comparison_AIC_BIC.csv")), row.names = FALSE)
  }
  
  
  #PCA & PERMANOVA
  cat("\n=== PCA & PERMANOVA ===\n")
  resp_mat <- analysis_data %>% dplyr::select(any_of(response_vars)) %>% dplyr::select(where(is.numeric)) %>% as.data.frame()
  if (ncol(resp_mat) >= 2) {
    pca_res <- prcomp(resp_mat, center = TRUE, scale. = TRUE)
    pca_plots <- create_pca_plots_comprehensive(pca_res, metadata = analysis_data, response_vars = response_vars)
    results$pca <- list(result=pca_res, plots=pca_plots)
    if (create_plots) { print(pca_plots$main); print(pca_plots$variance); print(pca_plots$loadings) }
    if (save_results) {
      ggsave(file.path("ordination", paste0(output_prefix, "_PCA_main.png")),    pca_plots$main,    width=9, height=7, dpi=300, bg="white")
      ggsave(file.path("ordination", paste0(output_prefix, "_PCA_variance.png")), pca_plots$variance, width=9, height=6, dpi=300, bg="white")
      ggsave(file.path("ordination", paste0(output_prefix, "_PCA_loadings.png")), pca_plots$loadings, width=9, height=7, dpi=300, bg="white")
      saveRDS(pca_res, file.path("ordination", paste0(output_prefix, "_PCA_result.rds")))
    }
  } else cat("PCA skipped (need ≥2 numeric response variables).\n")
  
  if (ncol(resp_mat) >= 1 && "Fishing_Ground" %in% names(analysis_data)) {
    dist_matrix <- vegan::vegdist(resp_mat, method = "euclidean")
    permanova_overall <- vegan::adonis2(dist_matrix ~ Fishing_Ground, data = analysis_data, permutations = 999)
    cat("\nOverall PERMANOVA by Fishing_Ground:\n"); print(permanova_overall)
    pairwise_perm <- pairwise.adonis2(resp_mat, analysis_data, permutations = 999, method = "euclidean", p_adjust = "BH")
    cat("\nPairwise PERMANOVA (BH):\n"); print(pairwise_perm$pairwise_results)
    results$permanova <- list(overall=permanova_overall, pairwise=pairwise_perm)
    if (save_results) {
      sink(file.path("tables", paste0(output_prefix, "_permanova_overall.txt"))); print(permanova_overall); sink()
      if (nrow(pairwise_perm$pairwise_results) > 0)
        write.csv(pairwise_perm$pairwise_results, file.path("tables", paste0(output_prefix, "_permanova_pairwise.csv")), row.names = FALSE)
    }
  } else cat("PERMANOVA skipped: insufficient response matrix or missing Fishing_Ground.\n")
  
  #Save LM models
  if (save_results && length(lm_results)) {
    lapply(names(lm_results), function(nm)
      saveRDS(lm_results[[nm]]$model, file.path("models", paste0(output_prefix, "_lm_", nm, ".rds"))))
  }
  
  #Univariate panel
  if (save_results) {
    univ_prefix <- paste0(output_prefix, "_univariate_", ifelse(run_univariate_scope=="all","ALL","excluded"))
    run_univariate_panel(
      analysis_data          = analysis_data,
      response_vars          = response_vars,
      all_predictors         = setdiff(unique(predictor_vars), "Fishing_Ground"),
      final_fixed_predictors = fixed_screened,
      random_effects         = random_effects,
      output_prefix          = univ_prefix,
      scope                  = run_univariate_scope,
      save_results           = TRUE
    )
  }
  
  #Partial Plots for Top Models
  if (create_conditional_partial_plots && save_results) {  
    cat("\n=== CREATING CONDITIONAL PARTIAL PLOTS FOR TOP MODELS ===\n")
    tryCatch({
      partial_plots <- create_conditional_partial_plots(results, output_prefix)  
      results$partial_plots <- partial_plots
      cat("Conditional partial plots completed!\n")
    }, error = function(e) {
      cat("Partial plots failed:", conditionMessage(e), "\n")
    })
  }
  
  
  #Report stub
  report_stub <- list(
    note = "Report generation placeholder for integration with R Markdown",
    created = Sys.time(),
    outputs = list(
      figures     = list.files("figures", full.names = TRUE),
      tables      = list.files("tables", full.names = TRUE),
      models      = list.files("models", full.names = TRUE),
      ordination  = list.files("ordination", full.names = TRUE),
      diagnostics = list.files("diagnostics", full.names = TRUE)
    ),
    screening = list(correlation_threshold=cor_threshold, vif_threshold=vif_threshold,
                     final_fixed_predictors=fixed_screened, scaled_predictors=scaled_cols)
  )
  results$report <- report_stub
  
  cat("\n", strrep("=", 74), "\n")
  cat("ANALYSIS COMPLETE. RESULTS OBJECT RETURNED.\n")
  cat(strrep("=", 74), "\n")
  return(results)
}


# EXECUTION OF ANALYSES


#Core ecological
cfg <- analysis_config$core_ecological
# prefix("core_ecological")
res_core <- analyze_reef_data_enhanced(
  response_vars = cfg$response_vars,
  predictor_vars = cfg$predictor_vars,
  random_effects = cfg$random_effects,
  data = clean_data,
  output_prefix = "core_ecological",
  create_plots = TRUE,
  save_results = TRUE,
  run_diagnostics = TRUE,
  run_nonparametric = TRUE,
  lmm_location_as_random_only = TRUE,
  scale_predictors_for_LMM = TRUE,
  scale_predictors_for_LM  = TRUE,
  cor_threshold = 0.75,
  vif_threshold = 5,
  run_univariate_scope = "excluded_only",
  create_conditional_partial_plots = TRUE  
)

#Functional diversity
cfg <- analysis_config$functional_diversity
# prefix("functionalDiv_traitbased")
res_fd <- analyze_reef_data_enhanced(
  response_vars = cfg$response_vars,
  predictor_vars = cfg$predictor_vars,
  random_effects = cfg$random_effects,
  data = clean_data,
  output_prefix = "functionalDiv_traitbased",
  create_plots = TRUE, save_results = TRUE, run_diagnostics = TRUE, run_nonparametric = TRUE,
  lmm_location_as_random_only = TRUE, scale_predictors_for_LMM = TRUE, scale_predictors_for_LM = TRUE,
  cor_threshold = 0.75, vif_threshold = 5, run_univariate_scope = "excluded_only",
  create_conditional_partial_plots = TRUE 
)

#Biomass focused
cfg <- analysis_config$biomass_focused
# purge_outputs_by_prefix("biomass_focused")
res_biomass <- analyze_reef_data_enhanced(
  response_vars = cfg$response_vars,
  predictor_vars = cfg$predictor_vars,
  random_effects = cfg$random_effects,
  data = clean_data,
  output_prefix = "biomass_focused",
  create_plots = TRUE, save_results = TRUE, run_diagnostics = TRUE, run_nonparametric = TRUE,
  lmm_location_as_random_only = TRUE, scale_predictors_for_LMM = TRUE, scale_predictors_for_LM = TRUE,
  cor_threshold = 0.75, vif_threshold = 5, run_univariate_scope = "excluded_only",
  create_conditional_partial_plots = TRUE 
)

#All responses
cfg <- analysis_config$all_responses
# purge_outputs_by_prefix("all_responses")
res_all <- analyze_reef_data_enhanced(
  response_vars = cfg$response_vars,
  predictor_vars = cfg$predictor_vars,
  random_effects = cfg$random_effects,
  data = clean_data,
  output_prefix = "all_responses",
  create_plots = TRUE, save_results = TRUE, run_diagnostics = TRUE, run_nonparametric = TRUE,
  lmm_location_as_random_only = TRUE, scale_predictors_for_LMM = TRUE, scale_predictors_for_LM = TRUE,
  cor_threshold = 0.75, vif_threshold = 5, run_univariate_scope = "excluded_only",
  create_conditional_partial_plots = TRUE  
)

