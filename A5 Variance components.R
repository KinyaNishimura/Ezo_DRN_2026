# Install and load required packages
packages <- c("RRPP", "dplyr", "lme4", "MuMIn", "ggplot2", "gridExtra", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Data loading
setwd('/Users/kin/Library/CloudStorage/GoogleDrive-kinya@eis.hokudai.ac.jp/マイドライブ/地域集団比較 2015 iCloud 211128 221027復帰/Anal Proc/TPS_and_other_files/５地域集団 解析/ReAnal 20240305 〜 4回/Anal BE')
data <- read.delim("Symm.BE.GPA.AllInd.txt", head=TRUE)

print("Data structure:")
print(str(data))
print("\nSample sizes by Location and Treatment:")
print(table(data$Location, data$Treatments))

##########     Step 1: PCA on GPA data
# GPA data processing
symm_cols <- paste0("Symm", 1:66)
gpa_data <- as.matrix(data[, symm_cols])

# Data centered around the mean
centered_data <- scale(gpa_data, center = TRUE, scale = FALSE)

# Execute PCA
pca_result <- prcomp(centered_data, scale. = FALSE)
print("\nPCA summary:")
print(summary(pca_result))

pc_scores <- pca_result$x
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cumulative_variance <- cumsum(var_explained)

print("\nPC variance proportions:")
pc_summary <- data.frame(
  PC = paste0("PC", 1:length(var_explained)),
  Variance_explained = var_explained,
  Cumulative = cumulative_variance
)
print(head(pc_summary, 15))

major_pcs <- which(var_explained >= 0.05 | cumulative_variance <= 0.95)
if(length(major_pcs) == 0) major_pcs <- 1:min(10, length(var_explained))

max_pc_to_analyze <- min(which(cumulative_variance >= 0.99), 20, length(var_explained))
cat("PCs to analyze:", max_pc_to_analyze, "\nMajor PCs:", paste(major_pcs, collapse = ", "), "\n")

##########     Step 2: Variance component analysis for each PC
pc_effects_analysis <- data.frame()

for (i in 1:max_pc_to_analyze) {
  pc_name <- paste0("PC", i)
  cat("Analyzing:", pc_name, "\n")
  
  pc_data <- data.frame(
    y = pc_scores[, i],
    Treatments = factor(data$Treatments),
    Location = factor(data$Location)
  )
  
  treatments_pvalue <- location_pvalue <- interaction_pvalue <- overall_pvalue <- NA
  treatments_ratio <- location_ratio <- interaction_ratio <- residual_ratio <- NA
  r2_marginal <- r2_conditional <- NA
  
  tryCatch({
    y_var <- var(pc_data$y)
    if(y_var < 1e-10) {
      cat("  Scaling: very small variance detected.\n")
      pc_data$y <- pc_data$y * 1e6
    }
    
    null_model <- lmer(y ~ 1 + (1|Location), data = pc_data, REML = FALSE,
                       control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
    treatments_model <- lmer(y ~ Treatments + (1|Location), data = pc_data, REML = FALSE,
                            control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
    
    full_model <- tryCatch({
      lmer(y ~ Treatments + (1|Location) + (1|Treatments:Location), data = pc_data, REML = FALSE,
           control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
    }, warning = function(w) {
      if(grepl("singular", w$message)) {
        cat("  Singular fit detected. Using simplified model.\n")
        lmer(y ~ Treatments + Location + (1|Location), data = pc_data, REML = FALSE,
             control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
      } else {
        lmer(y ~ Treatments + (1|Location) + (1|Treatments:Location), data = pc_data, REML = FALSE,
             control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
      }
    }, error = function(e) {
      cat("  Using linear model.\n")
      lm(y ~ Treatments + Location, data = pc_data)
    })
    
    is_lmer_model <- inherits(full_model, "lmerMod")
    
    if(is_lmer_model) {
      treatments_test <- tryCatch(anova(null_model, treatments_model), error = function(e) NULL)
      if(!is.null(treatments_test) && nrow(treatments_test) > 1) 
        treatments_pvalue <- treatments_test$`Pr(>Chisq)`[2]
      
      interaction_test <- tryCatch(anova(treatments_model, full_model), error = function(e) NULL)
      if(!is.null(interaction_test) && nrow(interaction_test) > 1) 
        interaction_pvalue <- interaction_test$`Pr(>Chisq)`[2]
      
      overall_test <- tryCatch(anova(null_model, full_model), error = function(e) NULL)
      if(!is.null(overall_test) && nrow(overall_test) > 1) 
        overall_pvalue <- overall_test$`Pr(>Chisq)`[nrow(overall_test)]
      
      vc <- tryCatch(as.data.frame(VarCorr(full_model)), error = function(e) data.frame(grp = character(), vcov = numeric()))
      
      location_var <- ifelse(any(vc$grp == "Location"), vc$vcov[vc$grp == "Location"][1], 0)
      interaction_var <- ifelse(any(vc$grp == "Treatments:Location"), vc$vcov[vc$grp == "Treatments:Location"][1], 0)
      residual_var <- ifelse(any(vc$grp == "Residual"), vc$vcov[vc$grp == "Residual"][1], 0)
      
      fitted_fixed <- tryCatch(predict(full_model, re.form = ~0), error = function(e) rep(0, nrow(pc_data)))
      treatments_var <- var(fitted_fixed)
      if(is.na(treatments_var)) treatments_var <- 0
      
      r2_results <- tryCatch(r.squaredGLMM(full_model), error = function(e) c(0, 0))
      r2_marginal <- r2_results[1]
      r2_conditional <- r2_results[2]
      
      location_model <- lm(y ~ Treatments, data = pc_data)
      location_test <- tryCatch(anova(location_model, lm(y ~ Treatments + Location, data = pc_data)), error = function(e) NULL)
      if(!is.null(location_test) && nrow(location_test) > 1) 
        location_pvalue <- location_test$`Pr(>F)`[2]
      
    } else {
      anova_result <- anova(full_model)
      if(nrow(anova_result) >= 2) {
        treatments_var <- anova_result$`Sum Sq`[1] / sum(anova_result$`Sum Sq`)
        location_var <- anova_result$`Sum Sq`[2] / sum(anova_result$`Sum Sq`)
        residual_var <- anova_result$`Sum Sq`[nrow(anova_result)] / sum(anova_result$`Sum Sq`)
        interaction_var <- 0
      }
      r2_marginal <- r2_conditional <- summary(full_model)$r.squared
    }
    
    total_var <- treatments_var + location_var + interaction_var + residual_var
    treatments_ratio <- ifelse(total_var > 0, treatments_var / total_var, 0)
    location_ratio <- ifelse(total_var > 0, location_var / total_var, 0)
    interaction_ratio <- ifelse(total_var > 0, interaction_var / total_var, 0)
    residual_ratio <- ifelse(total_var > 0, residual_var / total_var, 0)
    
    if(i %in% major_pcs) {
      cat(sprintf("  %s - PC variance components:\n", pc_name))
      cat(sprintf("    Treatments: %.1f%% | Location: %.1f%% | Interaction: %.1f%% | Residual: %.1f%%\n",
                  treatments_ratio*100, location_ratio*100, interaction_ratio*100, residual_ratio*100))
    }
    
  }, error = function(e) {
    cat("Error in", pc_name, ":", e$message, "\n")
  })
  
  pc_effects_analysis <- rbind(pc_effects_analysis, data.frame(
    PC = i,
    PC_name = pc_name,
    Variance_explained = var_explained[i],
    Cumulative_variance = cumulative_variance[i],
    Overall_pvalue = overall_pvalue,
    Treatments_pvalue = treatments_pvalue,
    Location_pvalue = location_pvalue,
    Interaction_pvalue = interaction_pvalue,
    Overall_significant = !is.na(overall_pvalue) && overall_pvalue < 0.05,
    Treatments_significant = !is.na(treatments_pvalue) && treatments_pvalue < 0.05,
    Location_significant = !is.na(location_pvalue) && location_pvalue < 0.05,
    Interaction_significant = !is.na(interaction_pvalue) && interaction_pvalue < 0.05,
    Treatments_ratio = treatments_ratio,
    Location_ratio = location_ratio,
    Interaction_ratio = interaction_ratio,
    Residual_ratio = residual_ratio,
    R2_marginal = r2_marginal,
    R2_conditional = r2_conditional,
    Weighted_importance = ifelse(!is.na(r2_conditional), r2_conditional * var_explained[i], 0)
  ))
}

print("\nPC variance component analysis completed")

##########     Step 3: Calculate contribution to total variance
effects_importance <- pc_effects_analysis %>%
  filter(PC %in% major_pcs) %>%
  mutate(
    Treatments_importance = Treatments_ratio * Variance_explained,
    Location_importance = Location_ratio * Variance_explained,
    Interaction_importance = Interaction_ratio * Variance_explained,
    Residual_importance = Residual_ratio * Variance_explained
  )

case1_results <- effects_importance %>%
  mutate(
    PC_total = Variance_explained * 100,
    Treatment = Treatments_importance * 100,
    Location = Location_importance * 100,
    Interaction = Interaction_importance * 100,
    Residual = Residual_importance * 100
  ) %>%
  select(PC = PC_name, PC_total, Treatment, Location, Interaction, Residual)

print("\nCase 1: Total variance = 100%")
print(case1_results)

non_residual_total <- sum(case1_results$Treatment + case1_results$Location + case1_results$Interaction)
case2_results <- case1_results %>%
  mutate(
    PC_total = ((Treatment + Location + Interaction) / non_residual_total) * 100,
    Treatment = (Treatment / non_residual_total) * 100,
    Location = (Location / non_residual_total) * 100,
    Interaction = (Interaction / non_residual_total) * 100
  ) %>%
  select(PC, PC_total, Treatment, Location, Interaction)

print("\nCase 2: Non-residual variance = 100%")
print(case2_results)

##########     Step 4: Visualizations
plot_data_case1 <- case1_results %>%
  pivot_longer(cols = c(Treatment, Location, Interaction, Residual), 
               names_to = "Component", values_to = "Percentage") %>%
  mutate(Component = factor(Component, levels = c("Treatment", "Location", "Interaction", "Residual")))

plot_data_case2 <- case2_results %>%
  pivot_longer(cols = c(Treatment, Location, Interaction), 
               names_to = "Component", values_to = "Percentage") %>%
  mutate(Component = factor(Component, levels = c("Treatment", "Location", "Interaction")))

colors_case1 <- c("Treatment" = "#66C2A5", "Location" = "#FC8D62", "Interaction" = "#8DA0CB", "Residual" = "#E78AC3")
colors_case2 <- c("Treatment" = "#66C2A5", "Location" = "#FC8D62", "Interaction" = "#8DA0CB")

p1 <- ggplot(plot_data_case1, aes(x = factor(PC, levels = paste0("PC", major_pcs)), y = Percentage, fill = Component)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors_case1) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Case 1: Total variance = 100%", x = "Principal Components", y = "Percentage (%)") +
  theme_minimal() + theme(legend.position = "top") + coord_flip()

p2 <- ggplot(plot_data_case2, aes(x = factor(PC, levels = paste0("PC", major_pcs)), y = Percentage, fill = Component)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors_case2) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Case 2: Non-residual variance = 100%", x = "Principal Components", y = "Percentage (%)") +
  theme_minimal() + theme(legend.position = "top") + coord_flip()

grid.arrange(p1, p2, ncol = 1)

cat("\nAnalysis completed.\n")
cat("Results available:\n")
cat("  - case1_results\n")
cat("  - case2_results\n")
cat("  - pc_effects_analysis\n")