## To create a correlation plot for the variables in TwinGene
library(BioAge)
# Define age variables
agevar <- c("PRS_Z", "age", "phenoage_res", "kdm_res", "hd_res", "Frailty_Index_0.10_log", "Telomere")
# Define axis types
axis_type <- c(
"PRS_Z" = "float",
"age" = "float",
"phenoage_res" = "float",
"kdm_res" = "float",
"hd_res" = "float",
"Frailty_Index_0.10_log" = "float",
"Telomere" = "float"
)
# Define labels with the specified renaming and order
label <- c(
"PRS_Z" = "Polygenic Risk Scores",
"age" = "Choronological age",
"phenoage_res" = "Phenoage\nLevine\nPhenotypic Age\nAdvancement",
"kdm_res" = "KDM\nModified-KDM\nBiological Age\nAdvancement",
"hd_res" = "HD\nHomeostatic\nDysregulation",
"Frailty_Index_0.10_log" = "Frailty Index",
"Telomere" = "Telomeres Length"
)
plot_baa(data_twingene, agevar, label, axis_type)
#---------------------------------------------------------------------------
## To calculate the univariate AUC (95% CI) for each predictor in the TwinGene
# Load necessary libraries
library(pROC)
library(dplyr)
# Ensure 'status' is a binary factor (0 = no event, 1 = event)
data_twingene$status <- as.factor(data_twingene$status)
# Define predictors
predictors <- c("age", "PRS_Z", "phenoage_res", "kdm_res", "hd_res",
"Frailty_Index_0.10_log", "Telomere")
# Function to compute AUC & 95% CI
calculate_auc <- function(var, data) {
roc_obj <- roc(data$status, data[[var]], ci = TRUE)
auc_value <- auc(roc_obj)
ci_values <- ci.auc(roc_obj)
return(data.frame(
Predictor = var,
AUC = round(auc_value, 3),
CI_Lower = round(ci_values[1], 3),
CI_Upper = round(ci_values[3], 3)
))
}
# Compute AUC for all predictors
auc_results <- do.call(rbind, lapply(predictors, calculate_auc, data = data_twingene))
# Compare each predictor with CA using DeLong’s test
p_values <- sapply(predictors[-1], function(var) {
roc_ca <- roc(data_twingene$status, data_twingene$age)
roc_var <- roc(data_twingene$status, data_twingene[[var]])
p_value <- roc.test(roc_ca, roc_var, method = "delong")$p.value
return(p_value)
})
# Format P-values and combine results
auc_results$p_value <- c("-", round(p_values, 10))  # Add p-values
# Print the final results table
print(auc_results)
#----------------------------------------------------------------------------------------------
## To compare AUCs of BA measures with CA
library(pROC)
library(dplyr)
# Ensure 'death_status' is a binary factor (0 = no event, 1 = event)
data_twingene$status <- as.factor(data_twingene$status)  # TwinGene dataset
# Define predictors
predictors <- c("PRS_Z", "phenoage_res", "kdm_res", "hd_res",
"Frailty_Index_0.10_log", "Telomere")
# Function to compute AUC & 95% CI
calculate_auc <- function(var, data, status_col) {
# Ensure there are no missing values for the predictor and the outcome
valid_data <- data %>%
filter(!is.na(.data[[var]]) & !is.na(.data[[status_col]]))
if (nrow(valid_data) == 0) {
stop(paste("No valid data for predictor:", var))
}
roc_obj <- roc(valid_data[[status_col]], valid_data[[var]], ci = TRUE)
auc_value <- auc(roc_obj)
ci_values <- ci.auc(roc_obj)
return(data.frame(
Predictor = var,
AUC = round(auc_value, 3),
CI_Lower = round(ci_values[1], 3),
CI_Upper = round(ci_values[3], 3)
))
}
# Compute AUC for all predictors in data_twingene
auc_results_twingene <- do.call(rbind, lapply(c("age", predictors), calculate_auc, data = data_twingene, status_col = "status"))
# Function to perform DeLong's test against CA
calculate_p_value <- function(var, data, status_col) {
roc_ca <- roc(data[[status_col]], data$age)
roc_var <- roc(data[[status_col]], data[[var]])
p_value <- roc.test(roc_ca, roc_var, method = "delong")$p.value
return(p_value)
}
# Calculate p-values for each predictor compared to CA in data_twingene
p_values_twingene <- sapply(predictors, calculate_p_value, data = data_twingene, status_col = "status")
# Combine AUC results and p-values for data_twingene
auc_results_twingene$p_value <- c("-", round(p_values_twingene, 10))  # Add p-values
# Print the final results table for Overall_data
print("Overall_data AUC Results:")
print(auc_results_overall)
# Print the final results table for TwinGene
print(auc_results_twingene)
#-------------------------------------------------------------------------------------------
## To create a forest plot of Univariate AUC for both cohorts 
library(ggplot2)
library(dplyr)
library(tidyr)
# Create the data frame with updated AUC values
data <- data.frame(
Predictor = c("CA", "PhenoAge", "KDM", "HD", "FI", "Telomere length", "PRS"),
AUC_TwinGene = c(0.835, 0.874, 0.747, 0.592, 0.574, 0.581, 0.510),
CI_Lower_TwinGene = c(0.819, 0.866, 0.736, 0.580, 0.562, 0.569, 0.498),
CI_Upper_TwinGene = c(0.848, 0.883, 0.758, 0.605, 0.587, 0.593, 0.522),
AUC_UKB = c(0.708, 0.624, 0.600, 0.565, 0.587, 0.521, 0.502),
CI_Lower_UKB = c(0.703, 0.619, 0.595, 0.560, 0.582, 0.516, 0.498),
CI_Upper_UKB = c(0.712, 0.629, 0.605, 0.570, 0.592, 0.526, 0.505)
)
# Order the data by AUC_TwinGene from highest to lowest
data <- data %>%
arrange(desc(AUC_TwinGene)) %>%
mutate(Predictor = factor(Predictor, levels = Predictor))
# Reshape the data for plotting
data_long <- data %>%
pivot_longer(cols = starts_with("AUC"), names_to = "Cohort", values_to = "AUC") %>%
mutate(Cohort = ifelse(grepl("TwinGene", Cohort), "TwinGene", "UK Biobank")) %>%
left_join(data %>%
pivot_longer(cols = starts_with("CI_Lower"), names_to = "Cohort", values_to = "CI_Lower") %>%
mutate(Cohort = ifelse(grepl("TwinGene", Cohort), "TwinGene", "UK Biobank")),
by = c("Predictor", "Cohort")) %>%
left_join(data %>%
pivot_longer(cols = starts_with("CI_Upper"), names_to = "Cohort", values_to = "CI_Upper") %>%
mutate(Cohort = ifelse(grepl("TwinGene", Cohort), "TwinGene", "UK Biobank")),
by = c("Predictor", "Cohort"))
# Plot using ggplot2 with enhanced font sizes and better colors
ggplot(data_long, aes(x = AUC, y = Predictor, color = Cohort)) +
geom_point(position = position_dodge(width = 0.5), size = 5) +  # Increase point size
geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper),
position = position_dodge(width = 0.5), height = 0.2, size = 1) +  # Thicker error bars
scale_color_manual(values = c("TwinGene" = "#D55E00", "UK Biobank" = "#0072B2")) +  # Distinct colors
labs(x = "AUC (95% CI)", y = "Predictor", color = "Cohort") +
theme_minimal() +
theme(
legend.position = "top",
legend.text = element_text(size = 16),  # Larger legend text
legend.title = element_text(size = 18, face = "bold"),
plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
axis.title = element_text(size = 20, face = "bold"),
axis.text = element_text(size = 18),
panel.grid.major = element_line(color = "gray85"),  # Subtle grid lines
panel.grid.minor = element_blank()
) +
ggtitle("Forest Plot of AUC with 95% CI for TwinGene and UK Biobank")
#----------------------------------------------------------------------------------------------
## To develop multivariate predictive model (including age and PRS) for mortality in the Twingene using Stacking ML
library(SuperLearner)
library(pROC)
library(caret)
library(randomForest)
library(xgboost)
# Define predictors excluding KDM and HD from the full model
predictors_age <- data_twingene[, "age", drop = FALSE]
predictors_reduced <- data_twingene[, c("age", "PRS_Z")]
predictors_full <- data_twingene[, c("age", "PRS_Z", "phenoage_res", "Telomere", "Frailty_Index_0.10_log", "gender", "bmi", "smoker", "education_years")]
outcome <- data_twingene$status
# Convert outcome to a numeric vector for classification
outcome <- as.numeric(as.character(outcome))
# Split data into training (70%) and testing (30%) sets
set.seed(123) # For reproducibility
train_index <- createDataPartition(outcome, p = 0.7, list = FALSE)
train_predictors_age <- predictors_age[train_index, ]
train_predictors_reduced <- predictors_reduced[train_index, ]
train_predictors_full <- predictors_full[train_index, ]
train_outcome <- outcome[train_index]
test_predictors_age <- predictors_age[-train_index, ]
test_predictors_reduced <- predictors_reduced[-train_index, ]
test_predictors_full <- predictors_full[-train_index, ]
test_outcome <- outcome[-train_index]
# Define base learners
base_learners <- c("SL.glm", "SL.rpart", "SL.randomForest", "SL.xgboost")
# Fit SuperLearner model on training data
sl_model_age <- SuperLearner(Y = train_outcome, X = train_predictors_age, family = binomial(), SL.library = base_learners, method = "method.NNLS")
sl_model_reduced <- SuperLearner(Y = train_outcome, X = train_predictors_reduced, family = binomial(), SL.library = base_learners, method = "method.NNLS")
sl_model_full <- SuperLearner(Y = train_outcome, X = train_predictors_full, family = binomial(), SL.library = base_learners, method = "method.NNLS")
# Predictions on testing data
predictions_age <- predict(sl_model_age, newdata = test_predictors_age, type = "response")
predictions_reduced <- predict(sl_model_reduced, newdata = test_predictors_reduced, type = "response")
predictions_full <- predict(sl_model_full, newdata = test_predictors_full, type = "response")
# Compute ROC and AUC
roc_curve_age <- roc(test_outcome, predictions_age$pred)
roc_curve_reduced <- roc(test_outcome, predictions_reduced$pred)
roc_curve_full <- roc(test_outcome, predictions_full$pred)
# Plot ROC curves with different curve types
plot(roc_curve_age, main = "ROC Curves Comparison in TwinGene", col = "#D55E00", lwd = 2, lty = 1, print.auc = TRUE, print.auc.y = 0.4, legacy.axes = TRUE)
plot(roc_curve_reduced, add = TRUE, col = "#D55E00", lwd = 2, lty = 2, print.auc = TRUE, print.auc.y = 0.3, legacy.axes = TRUE)
plot(roc_curve_full, add = TRUE, col = "#D55E00", lwd = 2, lty = 4, print.auc = TRUE, print.auc.y = 0.2, legacy.axes = TRUE)
legend("bottomright", legend = c("CA Model", "CA + PRS Model", "CA+PRS+BAs Model"), col = "#D55E00", lwd = 2, lty = c(1, 2, 4))
# Perform DeLong's test to compare ROC curves
delong_test <- roc.test(roc_curve_reduced, roc_curve_full, method = "delong")
print(paste("DeLong's Test p-value:", delong_test$p.value))
#----------------------------------------------------------------------------------------------
## Survival analyses
library(survival)
library(survminer)
library(rms)
variables <- c("age", "PRS_Z", "kdm_res", "phenoage_res", "hd_res", "Telomere", "Frailty_Index_0.10_log")
# Create an empty list to store results
univariate_results <- list()
for (var in variables) {
model <- coxph(Surv(time_to_event, status) ~ get(var), data = data_twingene)
summary_model <- summary(model)
# Extract HR and 95% CI
HR <- summary_model$coefficients[,"exp(coef)"]
lower_CI <- summary_model$conf.int[,"lower .95"]
upper_CI <- summary_model$conf.int[,"upper .95"]
p_value <- summary_model$coefficients[,"Pr(>|z|)"]
# Store results
univariate_results[[var]] <- data.frame(
Variable = var,
HR = HR,
Lower_95_CI = lower_CI,
Upper_95_CI = upper_CI,
P_value = p_value
)
}
# Combine all results into a dataframe
univariate_results_df <- do.call(rbind, univariate_results)
# Print results
print(univariate_results_df)
multivariable_model <- coxph(Surv(time_to_event, status) ~
age + PRS_Z + kdm_res + phenoage_res +
hd_res + Telomere + Frailty_Index_0.10_log+ gender + bmi +
smoker + education_years, data = data_twingene)
# Print the summary of the model
summary(multivariable_model)
c_index <- summary(multivariable_model)$concordance[1]
print(paste("C-index:", c_index))
###
library(survival)
library(survminer)
library(rms)
variables <- c("age", "PRS_Z", "kdm_res", "phenoage_res", "hd_res", "Telomere", "Frailty_Index_0.10_log")
# Create an empty list to store results
univariate_results <- list()
for (var in variables) {
model <- coxph(Surv(time_to_event, status) ~ get(var), data = data_twingene)
summary_model <- summary(model)
# Extract HR and 95% CI
HR <- summary_model$coefficients[,"exp(coef)"]
lower_CI <- summary_model$conf.int[,"lower .95"]
upper_CI <- summary_model$conf.int[,"upper .95"]
p_value <- summary_model$coefficients[,"Pr(>|z|)"]
# Store results
univariate_results[[var]] <- data.frame(
Variable = var,
HR = HR,
Lower_95_CI = lower_CI,
Upper_95_CI = upper_CI,
P_value = p_value
)
}
# Combine all results into a dataframe
univariate_results_df <- do.call(rbind, univariate_results)
# Print results
print(univariate_results_df)
multivariable_model <- coxph(Surv(time_to_event, status) ~
age + PRS_Z + phenoage_res +
Telomere + Frailty_Index_0.10_log+ gender + bmi +
smoker + education_years, data = data_twingene)
# Print the summary of the model
summary(multivariable_model)
c_index <- summary(multivariable_model)$concordance[1]
print(paste("C-index:", c_index))
#-----------------------------------------------------------------------
## ROC Curves for the TwinGene dataset
library(pROC)
library(ggplot2)
library(dplyr)
library(cowplot)
# Load your dataset (replace with actual file path or use readr::read_csv if needed)
# Ensure binary outcome is numeric (0 = alive, 1 = dead)
data_twingene$status <- as.numeric(data_twingene$status)
# Build logistic regression models
model_ca <- glm(status ~ age, data = data_twingene, family = binomial)
model_ca_prs <- glm(status ~ age + PRS, data = data_twingene, family = binomial)
model_full <- glm(status ~ age + PRS + phenoage_res + Frailty_Index_0.10_log + Telomere +
gender + bmi + smoker + education_years, data = data_twingene, family = binomial)
# Get predicted probabilities
data_twingene$pred_ca <- predict(model_ca, type = "response")
data_twingene$pred_ca_prs <- predict(model_ca_prs, type = "response")
data_twingene$pred_full <- predict(model_full, type = "response")
# Compute ROC curves
roc_ca <- roc(data_twingene$status, data_twingene$pred_ca)
roc_ca_prs <- roc(data_twingene$status, data_twingene$pred_ca_prs)
roc_full <- roc(data_twingene$status, data_twingene$pred_full)
# Plot ROC curves with light-to-dark color theme
ggroc(list(
`CA Model` = roc_ca,
`CA + PRS Model` = roc_ca_prs,
`CA + PRS + BAs Model` = roc_full
), aes = "color", size = 1.5) +
scale_color_manual(values = c(
"CA Model" = "#FDBE85",           # Light orange
"CA + PRS Model" = "#FD8D3C",     # Medium orange
"CA + PRS + BAs Model" = "#E6550D"  # Dark orange
)) +
labs(title = "ROC Curves Comparison in TwinGene",
x = "1 - Specificity", y = "Sensitivity",
color = "") +
annotate("text", x = 0.6, y = 0.35, label = paste("AUC:", round(auc(roc_ca), 3)), color = "#FDBE85", size = 5) +
annotate("text", x = 0.6, y = 0.30, label = paste("AUC:", round(auc(roc_ca_prs), 3)), color = "#FD8D3C", size = 5) +
annotate("text", x = 0.6, y = 0.25, label = paste("AUC:", round(auc(roc_full), 3)), color = "#E6550D", size = 5) +
theme_cowplot(font_size = 14) +
theme(legend.position = c(0.75, 0.2),
plot.title = element_text(hjust = 0.5, face = "bold"))
#--------------------------------------------------------------------------------------------
## Forest plot of HRs in both TwinGene and UKB
library(ggplot2)
library(dplyr)
library(tidyr)
# Create the data frame
data <- data.frame(
Predictor = c("CA", "PhenoAge", "KDM", "HD", "FI", "Telomere length", "PRS"),
HR_TwinGene = c(1.15, 1.30, 1.16, 1.29, 1.33, 0.66, 0.99),
CI_Lower_TwinGene = c(1.14, 1.29, 1.15, 1.26, 1.26, 0.62, 0.96),
CI_Upper_TwinGene = c(1.16, 1.31, 1.17, 1.33, 1.43, 0.71, 1.03),
HR_UKB = c(1.102, 1.06, 1.068, 1.21, 1.57, 0.58, 0.996),
CI_Lower_UKB = c(1.10, 1.056, 1.065, 1.19, 1.55, 0.52, 0.990),
CI_Upper_UKB = c(1.105, 1.062, 1.071, 1.22, 1.60, 0.64, 1.010),
C_Index_TwinGene = 0.913,
C_Index_UKB = 0.749
)
# Order data
data <- data %>%
arrange(desc(HR_TwinGene)) %>%
mutate(Predictor = factor(Predictor, levels = Predictor))
# Reshape the data to long format
data_long <- data %>%
pivot_longer(cols = c(HR_TwinGene, HR_UKB), names_to = "Cohort", values_to = "HR") %>%
mutate(
Cohort = ifelse(grepl("TwinGene", Cohort), "TwinGene", "UK Biobank"),
CI_Lower = ifelse(Cohort == "TwinGene", CI_Lower_TwinGene, CI_Lower_UKB),
CI_Upper = ifelse(Cohort == "TwinGene", CI_Upper_TwinGene, CI_Upper_UKB)
)
# Define C-indices as subtitle
c_index_text <- paste0("C-index (TwinGene): ", data$C_Index_TwinGene[1],
" | C-index (UK Biobank): ", data$C_Index_UKB[1])
# Plot the data with increased font size
ggplot(data_long, aes(x = HR, y = Predictor, color = Cohort)) +
geom_point(position = position_dodge(width = 0.5), size = 3) +
geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, position = position_dodge(width = 0.5)) +
scale_color_manual(values = c("TwinGene" = "#D55E00", "UK Biobank" = "#0072B2")) +
labs(
x = "Hazard Ratio (95% CI)",
y = "Predictor",
color = "Cohort",
title = "Forest Plot of Hazard Ratios with 95% CI for TwinGene and UK Biobank",
subtitle = c_index_text
) +
theme_minimal() +
theme(
legend.position = "top",
text = element_text(size = 18),
plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
plot.subtitle = element_text(size = 16, hjust = 0.5, vjust = 2),
axis.title = element_text(size = 18),
axis.text = element_text(size = 16)
)
#------------------------------------------------------------------------------------------------------------------------
## Univariate AUCs for each predictor in the TwinGene with subgroup analysis for age categories and follow-up time
library(pROC)
library(dplyr)
# Ensure 'status' is a binary factor (0 = no event, 1 = event)
data_twingene$status <- as.factor(data_twingene$status)
# calculation for death_status in each follow-up period
data_twingene <- data_twingene %>%
mutate(
death_status_5y  = ifelse(time_to_event < 5  & status == 1, 1, ifelse(time_to_event >= 5, 0, NA)),
death_status_10y = ifelse(time_to_event < 10 & status == 1, 1, ifelse(time_to_event >= 10, 0, NA)),
death_status_15y = ifelse(time_to_event < 15 & status == 1, 1, ifelse(time_to_event >= 15, 0, NA))
)
# Categorize age into three groups
data_twingene <- data_twingene %>%
mutate(
age_group = case_when(
age < 50 ~ "<50",
age >= 50 & age <= 60 ~ "50-60",
age > 60 ~ ">60"
)
)
# Define predictors
predictors <- c("PRS_Z", "phenoage_res", "kdm_res", "hd_res",
"Frailty_Index_0.10_log", "Telomere")
# Function to compute AUC & 95% CI for subgroups
calculate_auc <- function(var, data, status_col) {
valid_data <- data %>%
filter(!is.na(.data[[var]]) & !is.na(.data[[status_col]]))
if (nrow(valid_data) == 0) {
return(data.frame(Predictor = var, AUC = NA, CI_Lower = NA, CI_Upper = NA))
}
roc_obj <- roc(valid_data[[status_col]], valid_data[[var]], ci = TRUE)
auc_value <- auc(roc_obj)
ci_values <- ci.auc(roc_obj)
return(data.frame(
Predictor = var,
AUC = round(auc_value, 3),
CI_Lower = round(ci_values[1], 3),
CI_Upper = round(ci_values[3], 3)
))
}
# Function to compute AUCs for subgroups
compute_auc_subgroups <- function(data, followup_col, age_group) {
subset_data <- data %>% filter(age_group == age_group)
auc_results <- do.call(rbind, lapply(predictors, calculate_auc, data = subset_data, status_col = followup_col))
auc_results$Age_Group <- age_group
auc_results$Follow_Up <- followup_col
return(auc_results)
}
# Perform subgroup analysis
subgroup_results <- bind_rows(
compute_auc_subgroups(data_twingene, "death_status_5y", "<50"),
compute_auc_subgroups(data_twingene, "death_status_5y", "50-60"),
compute_auc_subgroups(data_twingene, "death_status_5y", ">60"),
compute_auc_subgroups(data_twingene, "death_status_10y", "<50"),
compute_auc_subgroups(data_twingene, "death_status_10y", "50-60"),
compute_auc_subgroups(data_twingene, "death_status_10y", ">60"),
compute_auc_subgroups(data_twingene, "death_status_15y", "<50"),
compute_auc_subgroups(data_twingene, "death_status_15y", "50-60"),
compute_auc_subgroups(data_twingene, "death_status_15y", ">60")
)
# Print the final results table
print(subgroup_results)
#--------------------------------------------------------------------------------------------------------------------------
## To check interaction between BMI and BAs performance in TwinGene
library(pROC)
# Ensure 'status' is a binary factor (0 = no event, 1 = event)
data_twingene$status <- as.factor(data_twingene$status)
# Function to calculate AUC and p-value for interaction term between BMI and predictor
calculate_auc_interaction_with_pvalue <- function(predictor, data, status_col) {
# Remove rows with missing values for BMI and the predictor
valid_data <- data[!is.na(data$bmi) & !is.na(data[[predictor]]) & !is.na(data[[status_col]]), ]
# Create interaction term between BMI and the predictor
valid_data$interaction_term <- valid_data$bmi * valid_data[[predictor]]
# Fit logistic regression model with interaction term
formula <- as.formula(paste(status_col, "~ bmi +", predictor, "+ interaction_term"))
model <- glm(formula, data = valid_data, family = binomial)
# Get predicted probabilities from the model
predicted_probs <- predict(model, type = "response")
# Compute ROC curve
roc_obj <- roc(valid_data[[status_col]], predicted_probs, ci = TRUE)
# Get AUC and its confidence interval
auc_value <- auc(roc_obj)
ci_values <- ci.auc(roc_obj)
# Get p-value for the interaction term (bmi:predictor)
interaction_pvalue <- summary(model)$coefficients["interaction_term", "Pr(>|z|)"]
# Return the results with AUC, CI, and p-value for the interaction term
return(data.frame(
Predictor = predictor,
AUC = round(auc_value, 3),
CI_Lower = round(ci_values[1], 3),
CI_Upper = round(ci_values[3], 3),
Interaction_PValue = round(interaction_pvalue, 4)
))
}
# Define predictors
predictors <- c("PRS_Z", "phenoage_res", "kdm_res", "hd_res",
"Frailty_Index_0.10_log", "Telomere")
# Calculate AUCs and p-values for each predictor and its interaction with BMI
auc_results_interaction <- do.call(rbind, lapply(predictors, calculate_auc_interaction_with_pvalue, data = data_twingene, status_col = "status"))
# Print AUC results with interaction terms and p-values
print(auc_results_interaction)
#-----------------------------------------------------------------------------------------------------------------------------------------------
## Sub-group analysis of predictive power of each predictor between sexes in TwinGene
library(pROC)
data_twingene$status <- as.factor(data_twingene$status)
# Function to calculate AUC for a predictor, stratified by gender
calculate_auc_by_gender <- function(predictor, data) {
# Filter data by gender: 1 = Men, 2 = Women
data_men <- data[data$gender == 1, ]
data_women <- data[data$gender == 2, ]
# Remove NA values for the predictor
data_men <- data_men[!is.na(data_men[[predictor]]), ]
data_women <- data_women[!is.na(data_women[[predictor]]), ]
# Compute ROC for men and women
roc_men <- roc(data_men$status, data_men[[predictor]], ci = TRUE)
roc_women <- roc(data_women$status, data_women[[predictor]], ci = TRUE)
# Compute AUC values
auc_men <- auc(roc_men)
auc_women <- auc(roc_women)
# Perform DeLong's test to compare AUC between genders
test_result <- roc.test(roc_men, roc_women, method = "delong")
# Extract p-value from DeLong’s test
p_value <- test_result$p.value
# Return results as a dataframe
return(data.frame(
Predictor = predictor,
AUC_Men = round(auc_men, 3),
AUC_Women = round(auc_women, 3),
P_Value = round(p_value, 4)
))
}
# Define predictors to test
predictors <- c("PRS_Z", "phenoage_res", "kdm_res", "hd_res",
"Frailty_Index_0.10_log", "Telomere")
# Compute AUCs and statistical comparisons for each predictor
auc_gender_results <- do.call(rbind, lapply(predictors, calculate_auc_by_gender, data = data_twingene))
# Print the results
print(auc_gender_results)
## in the UKB
library(pROC)
# Ensure 'death_status' is a binary factor (0 = no event, 1 = event)
Overall_data$death_status <- as.factor(Overall_data$death_status)
# Function to calculate AUC for a predictor, stratified by gender
calculate_auc_by_gender <- function(predictor, data) {
# Filter data by gender: 1 = Men, 0 = Women
data_men <- data[data$sex == 1, ]   # Sex = 1 for Men
data_women <- data[data$sex == 0, ] # Sex = 0 for Women
# Remove NA values for the predictor
data_men <- data_men[!is.na(data_men[[predictor]]), ]
data_women <- data_women[!is.na(data_women[[predictor]]), ]
# Compute ROC for men and women
roc_men <- roc(data_men$death_status, data_men[[predictor]], ci = TRUE)
roc_women <- roc(data_women$death_status, data_women[[predictor]], ci = TRUE)
# Compute AUC values
auc_men <- auc(roc_men)
auc_women <- auc(roc_women)
# Perform DeLong's test to compare AUC between genders
test_result <- roc.test(roc_men, roc_women, method = "delong")
# Extract p-value from DeLong’s test
p_value <- test_result$p.value
# Return results as a dataframe
return(data.frame(
Predictor = predictor,
AUC_Men = round(auc_men, 3),
AUC_Women = round(auc_women, 3),
P_Value = round(p_value, 4)
))
}
# Define predictors to test (adjust these if needed)
predictors <- c("PRS_Z", "Phenoage.Residual", "HD.Residual", "KDM.Residual",
"FI.Residual", "Telomere_Length.Residual")
# Compute AUCs and statistical comparisons for each predictor
auc_gender_results <- do.call(rbind, lapply(predictors, calculate_auc_by_gender, data = Overall_data))
# Print the results
print(auc_gender_results)
#-------------------------------------------------------------------------------------------------------------------
## Sub-group analysis of predictive power of each predictor between smokers and non-smokers in TwinGene
library(pROC)
# Ensure 'status' is a binary factor (0 = no event, 1 = event)
data_twingene$status <- as.factor(data_twingene$status)
# Function to calculate AUC for a predictor, stratified by smoker status
calculate_auc_by_smoker <- function(predictor, data) {
# Filter data by smoker status: 0 = Non-Smoker, 1 = Smoker
data_non_smoker <- data[data$smoker == 0, ]
data_smoker <- data[data$smoker == 1, ]
# Remove NA values for the predictor
data_non_smoker <- data_non_smoker[!is.na(data_non_smoker[[predictor]]), ]
data_smoker <- data_smoker[!is.na(data_smoker[[predictor]]), ]
# Compute ROC for non-smokers and smokers
roc_non_smoker <- roc(data_non_smoker$status, data_non_smoker[[predictor]], ci = TRUE)
roc_smoker <- roc(data_smoker$status, data_smoker[[predictor]], ci = TRUE)
# Compute AUC values
auc_non_smoker <- auc(roc_non_smoker)
auc_smoker <- auc(roc_smoker)
# Perform DeLong's test to compare AUC between smoker and non-smoker groups
test_result <- roc.test(roc_non_smoker, roc_smoker, method = "delong")
# Extract p-value from DeLong’s test
p_value <- test_result$p.value
# Return results as a dataframe
return(data.frame(
Predictor = predictor,
AUC_Non_Smoker = round(auc_non_smoker, 3),
AUC_Smoker = round(auc_smoker, 3),
P_Value = round(p_value, 4)
))
}
# Define predictors to test
predictors <- c("PRS_Z", "phenoage_res", "kdm_res", "hd_res",
"Frailty_Index_0.10_log", "Telomere")
# Compute AUCs and statistical comparisons for each predictor
auc_smoker_results <- do.call(rbind, lapply(predictors, calculate_auc_by_smoker, data = data_twingene))
# Print the results
print(auc_smoker_results)
## in the UKB
library(pROC)
# Ensure 'death_status' is a binary factor (0 = no event, 1 = event)
Overall_data$death_status <- as.factor(Overall_data$death_status)
# Function to calculate AUC for a predictor, stratified by smoker status
calculate_auc_by_smoker <- function(predictor, data) {
# Filter data by smoker status: 1 = Non-Smoker, 2 = Smoker
data_non_smoker <- data[data$smoking == 1, ]  # Non-smokers
data_smoker <- data[data$smoking == 2, ]      # Smokers
# Remove NA values for the predictor
data_non_smoker <- data_non_smoker[!is.na(data_non_smoker[[predictor]]), ]
data_smoker <- data_smoker[!is.na(data_smoker[[predictor]]), ]
# Compute ROC for non-smokers and smokers
roc_non_smoker <- roc(data_non_smoker$death_status, data_non_smoker[[predictor]], ci = TRUE)
roc_smoker <- roc(data_smoker$death_status, data_smoker[[predictor]], ci = TRUE)
# Compute AUC values
auc_non_smoker <- auc(roc_non_smoker)
auc_smoker <- auc(roc_smoker)
# Perform DeLong's test to compare AUC between smoker and non-smoker groups
test_result <- roc.test(roc_non_smoker, roc_smoker, method = "delong")
# Extract p-value from DeLong’s test
p_value <- test_result$p.value
# Return results as a dataframe
return(data.frame(
Predictor = predictor,
AUC_Non_Smoker = round(auc_non_smoker, 3),
AUC_Smoker = round(auc_smoker, 3),
P_Value = round(p_value, 4)
))
}
# Define predictors to test (update as needed)
predictors <- c("PRS_Z", "Phenoage.Residual", "HD.Residual", "KDM.Residual",
"FI.Residual", "Telomere_Length.Residual")
# Compute AUCs and statistical comparisons for each predictor
auc_smoker_results <- do.call(rbind, lapply(predictors, calculate_auc_by_smoker, data = Overall_data))
# Print the results
print(auc_smoker_results)
#-------------------------------------------------------------------------------------------------------------------------------------
## Sub-group analysis of predictive power of BAs and education level in TwinGene
# Load necessary library
library(dplyr)
# Create the interaction model
model <- glm(status ~ education_years * PRS_Z + education_years * phenoage_res + education_years * hd_res + education_years * kdm_res +
education_years * Telomere + education_years * Frailty_Index_0.10_log,
family = binomial,
data = data_twingene)
# Display the summary of the model
summary(model)
#--------------------------------------------------------------------------------------------------------
## S.Figure 1
library(tidyverse)
# Step 1: Prepare transformed variables with new labels
data_prepped <- data_twingene %>%
transmute(
`glucose (log-transformed; mmol/L)` = log(glucose_mmol),
`high-density lipoprotein (HDL)` = hdl,
`low-density lipoprotein (LDL)` = ldl,
`total cholesterol (mg/dL)` = totchol * 38.67,
`triglycerides (log-transformed; mg/dL)` = log(trig * 88.57),
`hemoglobin A1c (HbA1c; log-transformed)` = lnhba1c,
`creatinine (log-transformed; μmol/L)` = lncreat_umol,
`cystatin C (mg/L; untransformed)` = cyst,
`C-reactive protein (CRP; log-transformed; mg/dL)` = log(exp(trig) / 10)  # placeholder
)
# Step 2: Reshape to long format
df_long <- data_prepped %>%
pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
filter(!is.na(Value))
# Step 3: Create data for overlaying normal curves within each facet
normal_curves <- df_long %>%
group_by(Variable) %>%
summarise(mean = mean(Value), sd = sd(Value), .groups = "drop") %>%
rowwise() %>%
mutate(x = list(seq(mean - 4 * sd, mean + 4 * sd, length.out = 200))) %>%
unnest(x) %>%
mutate(density = dnorm(x, mean, sd))
# Step 4: Final plot
ggplot(df_long, aes(x = Variable, y = Value)) +
geom_violin(fill = "#D55E00", color = "black", alpha = 0.6, trim = FALSE) +
stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +
geom_line(data = normal_curves,
aes(x = Variable, y = x, group = Variable),  # x = Variable for facetting
inherit.aes = FALSE, color = "blue", linetype = "dashed") +
facet_wrap(~Variable, scales = "free", ncol = 3) +
labs(title = "TwinGene Biomarker Distributions", x = "", y = "Value") +
theme_minimal(base_size = 14) +
theme(axis.text.x = element_blank(),
strip.text = element_text(face = "bold"),
legend.position = "none")
#------------------------------------------------------------------------------------------------------------------
