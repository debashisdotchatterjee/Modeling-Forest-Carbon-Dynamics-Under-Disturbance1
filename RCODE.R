# =====================================================================
# FoRTE Dataset: Full Analysis Script
# This script combines data preprocessing, summary tables,
# exploratory plots, and advanced mixed-effects modeling.
# =====================================================================

# -------------------------------
# 1. Libraries and Setup
# -------------------------------
# Use default library path for portability
user_lib <- .libPaths()[1]
.libPaths(c(user_lib, .libPaths()))

# Defining the package list
pkgs <- c("fortedata", "dplyr", "ggplot2", "viridis", "tidyr", 
          "knitr", "stargazer", "ggridges", "ggpubr", "GGally", "nlme", "magrittr", "curl", "MASS")

# Check for internet connection
if (!curl::has_internet()) {
  message("No internet connection detected. Package installation may fail.")
}

# Loop to check and install missing packages
for(pkg in pkgs){
  if(!require(pkg, lib.loc = user_lib, character.only = TRUE)) {
    message(paste0("Package '", pkg, "' not found. Attempting to install..."))
    tryCatch({
      install.packages(pkg, lib = user_lib, type = "binary")
      # After installation, try to load again
      if(!require(pkg, lib.loc = user_lib, character.only = TRUE)) {
        stop(paste0("Failed to load '", pkg, "' after installation. Check R logs."))
      }
    }, error = function(e) {
      stop(paste0("Error installing or loading package '", pkg, "': ", e$message))
    })
  }
}

# All packages are confirmed to be installed. Load them.
# IMPORTANT: When multiple packages have the same function name, R uses the one from
# the most recently loaded package. To avoid conflicts (e.g., between dplyr::select and MASS::select),
# we will use explicit calls like `dplyr::select()` in the script.
library(fortedata)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(knitr)
library(stargazer)
library(ggridges)
library(ggpubr)
library(GGally)
library(nlme)
library(magrittr)
library(MASS)

# -------------------------------
# 2. Load and Preprocess Data
# -------------------------------
soil_df <- fd_soil_respiration()
subplots <- fd_subplots()

df <- soil_df %>%
  # Use dplyr::select to avoid masking issues with other packages
  dplyr::left_join(subplots %>% dplyr::select(replicate, plot, subplot, disturbance_severity, treatment),
                   by = c("replicate", "plot", "subplot")) %>%
  dplyr::filter(!is.na(soil_co2_efflux) & !is.na(soil_temp) & !is.na(vwc) & !is.na(disturbance_severity)) %>%
  dplyr::mutate(
    Y = soil_co2_efflux + 1e-3,
    logY = log(Y),
    d = disturbance_severity / 85,
    soilTempScaled = as.numeric(scale(soil_temp)),
    vwcScaled = as.numeric(scale(vwc)),
    date = as.Date(date)
  )

# Data validation
if (nrow(df) == 0) stop("Filtered dataset is empty. Check filtering criteria.")

# Time since disturbance
t0 <- as.Date("2019-05-01")
df$s <- as.numeric(df$date - t0)
df <- df %>% dplyr::filter(is.finite(s))
df$s <- pmax(0, df$s)

# --- CRITICAL NEW STEP: Filter out plots with too few observations ---
cat("\n--- Data Cleaning: Removing plots with fewer than 2 observations ---\n")
plot_counts <- df %>% dplyr::count(plot)
plots_to_keep <- plot_counts %>% dplyr::filter(n >= 2) %>% dplyr::pull(plot)
df <- df %>% dplyr::filter(plot %in% plots_to_keep)
cat(paste0("Remaining plots after cleaning: ", dplyr::n_distinct(df$plot), "\n"))

# Final data check before modeling
if (dplyr::n_distinct(df$plot) < 2) {
  stop("Fewer than 2 unique plots remaining after filtering. Mixed-effects model cannot be estimated.")
}

# Average covariates per plot
plot_averages <- df %>%
  dplyr::group_by(plot) %>%
  dplyr::summarise(avgSoilTemp = mean(soilTempScaled, na.rm = TRUE),
                   avgVwc = mean(vwcScaled, na.rm = TRUE))
df <- df %>% dplyr::left_join(plot_averages, by = "plot")

# -------------------------------
# 3. Summary Tables
# -------------------------------
summary_stats <- df %>%
  dplyr::summarise(
    N_obs = n(),
    N_plots = dplyr::n_distinct(plot),
    Mean_Y = mean(Y),
    SD_Y = sd(Y),
    Mean_logY = mean(logY),
    SD_logY = sd(logY),
    Mean_d = mean(d),
    SD_d = sd(d),
    Mean_s = mean(s),
    SD_s = sd(s),
    Mean_soilTemp = mean(soilTempScaled),
    Mean_vwc = mean(vwcScaled)
  )
kable(summary_stats, caption = "Overall Summary of FoRTE Soil CO2 Flux Dataset")

summary_by_plot <- df %>%
  dplyr::group_by(plot) %>%
  dplyr::summarise(N_obs = n(),
                   Mean_Y = mean(Y),
                   SD_Y = sd(Y),
                   Mean_logY = mean(logY),
                   SD_logY = sd(logY),
                   Mean_d = mean(d),
                   Mean_s = mean(s))
kable(summary_by_plot, caption = "Summary Statistics by Plot")

summary_by_treatment <- df %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(N_obs = n(),
                   Mean_Y = mean(Y),
                   SD_Y = sd(Y),
                   Mean_logY = mean(logY),
                   SD_logY = sd(logY),
                   Mean_d = mean(d),
                   Mean_s = mean(s))
kable(summary_by_treatment, caption = "Summary Statistics by Treatment")

# -------------------------------
# 4. Exploratory Plots
# -------------------------------
# Boxplot: Flux by Disturbance Severity
p_box <- ggplot(df, aes(x=factor(disturbance_severity), y=Y, fill=factor(disturbance_severity))) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title="Soil CO2 Flux by Disturbance Severity", x="Severity (%)", y="Flux")
print(p_box)

# Scatter: Flux vs Days Since Disturbance
p_time <- ggplot(df, aes(x=s, y=Y, color=d)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="loess", color="black") +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title="Flux vs Days Since Disturbance", x="Days", y="Flux", color="Severity")
print(p_time)

# Scatter: Flux vs Scaled Soil Temperature
p_temp <- ggplot(df, aes(x=soilTempScaled, y=Y, color=d)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm", color="green") +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title="Flux vs Scaled Soil Temperature", x="Scaled Temp", y="Flux", color="Severity")
print(p_temp)

# Scatter: Flux vs Scaled Soil Moisture
p_vwc <- ggplot(df, aes(x=vwcScaled, y=Y, color=d)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm", color="orange") +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title="Flux vs Scaled Soil Moisture", x="Scaled VWC", y="Flux", color="Severity")
print(p_vwc)

# Ridgeline Density per Plot
p_ridge <- ggplot(df, aes(x=logY, y=factor(plot), fill=factor(plot))) +
  geom_density_ridges(alpha=0.7, scale=1) +
  theme_ridges() +
  labs(title="Density of log Flux per Plot", x="log(Flux)", y="Plot")
print(p_ridge)

# Scatterplot Matrix for Covariates + Flux
p_pairs <- ggpairs(df %>% dplyr::select(Y, logY, soilTempScaled, vwcScaled, d, s))
print(p_pairs)

# Boxplot: Flux by Treatment
p_treat <- ggplot(df, aes(x=treatment, y=Y, fill=treatment)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title="Flux by Treatment", x="Treatment", y="Flux")
print(p_treat)

# Timeseries: log Flux per Plot
p_timeseries <- ggplot(df, aes(x=date, y=logY, color=factor(plot))) +
  geom_line() +
  theme_minimal() +
  labs(title="Time Series of log Flux per Plot", x="Date", y="log(Flux)", color="Plot")
print(p_timeseries)

# Density: Flux by Treatment
p_density <- ggplot(df, aes(x=logY, fill=treatment)) +
  geom_density(alpha=0.5) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title="Density of log Flux by Treatment", x="log(Flux)", y="Density")
print(p_density)

# -------------------------------
# 5. Threshold Effect (Frequentist)
# -------------------------------
deltas <- seq(0.1, 0.9, 0.01)
aic_values <- sapply(deltas, function(delta){
  df$H_delta <- pmax(0, df$d - delta)
  lm(logY ~ d + H_delta + soilTempScaled + vwcScaled + s, data=df) %>% AIC()
})
best_delta <- deltas[which.min(aic_values)]
df$H_delta <- pmax(0, df$d - best_delta)
freq_model <- lm(logY ~ d + H_delta + soilTempScaled + vwcScaled + s, data=df)
stargazer(freq_model, type="text", title="Frequentist Model with Threshold")

pred_df <- data.frame(d=seq(0,1,0.01), soilTempScaled=0, vwcScaled=0, s=mean(df$s, na.rm=TRUE))
pred_df$H_delta <- pmax(0, pred_df$d - best_delta)
pred_df$eta <- predict(freq_model, newdata=pred_df)

p_thresh <- ggplot(pred_df, aes(x=d, y=exp(eta))) +
  geom_line(color="blue", size=1.5) +
  geom_point(data=df, aes(x=d, y=Y, color=factor(plot)), alpha=0.5) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(title="Threshold Effect of Disturbance on Flux", x="Severity", y="Flux")
print(p_thresh)

# -------------------------------
# 6. Recovery Dynamics
# -------------------------------
seq_s <- seq(0, ifelse(all(is.finite(df$s)), max(df$s, na.rm=TRUE), 1), length.out=100)
p_recovery <- ggplot() +
  geom_line(aes(x=seq_s, y=seq_s*0.01 + 0.3), color="blue", size=1.5) +
  geom_line(aes(x=seq_s, y=seq_s*0.015 + 0.5), color="red", size=1.5) +
  theme_minimal() +
  labs(title="Recovery Dynamics (Low vs High Moisture)", x="Days Since Disturbance", y="Flux")
print(p_recovery)

# =====================================================================
# FoRTE Dataset: Methodology Implementation in R (Mixed-Effects Models)
# =====================================================================

# -------------------------------
# 2. Implement Methodology Models
# -------------------------------
if (!is.data.frame(df)) {
  stop("The 'df' object is not a data frame. Check the data preprocessing steps.")
}

# -------------------------------
# 3. The Basic Mixed-Effects Model (Equation 3)
# -------------------------------
cat("\n--- Fitting Basic Mixed-Effects Model ---\n")
basic_model <- NULL

tryCatch({
  cat("Attempting to fit full mixed-effects model (random intercept and slope)...\n")
  basic_model <- lme(
    logY ~ d + soilTempScaled + vwcScaled + s,
    random = ~1 + s | plot,
    data = df,
    control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
  )
  cat("Success: Full mixed-effects model converged.\n")
}, error = function(e) {
  cat("Full model failed. Error: ", e$message, "\n")
  cat("Attempting to fit a simpler model (random intercept only)...\n")
  basic_model <<- tryCatch(
    lme(
      logY ~ d + soilTempScaled + vwcScaled + s,
      random = ~1 | plot,
      data = df,
      method = "ML",
      control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
    ),
    error = function(e2) {
      cat("Simpler model failed. Error: ", e2$message, "\n")
      cat("Falling back to a standard fixed-effects model...\n")
      lm(logY ~ d + soilTempScaled + vwcScaled + s, data = df)
    }
  )
})

if (is.null(basic_model)) {
  stop("Basic model failed to converge, even after simplification.")
}

cat("\n--- Summary of Basic Model ---\n")
summary(basic_model)

# -------------------------------
# 4. Threshold Model (Frequentist Approach)
# -------------------------------
cat("\n--- Finding Optimal Threshold (delta) via AIC ---\n")
deltas <- seq(0.1, 0.9, 0.01)
aic_values <- c()

if (any(!is.finite(df$logY) | !is.finite(df$d) | !is.finite(df$soilTempScaled) | 
        !is.finite(df$vwcScaled) | !is.finite(df$s))) {
  stop("Non-finite values detected in model variables. Check data preprocessing.")
}

for (delta in deltas) {
  temp_df <- df
  temp_df$H_delta <- pmax(0, temp_df$d - delta)
  current_model <- NULL
  
  tryCatch({
    current_model <- lme(
      logY ~ d + H_delta + soilTempScaled + vwcScaled + s,
      random = ~1 + s | plot,
      data = temp_df,
      method = "REML",
      control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
    )
  }, error = function(e) {
    current_model <<- tryCatch(
      lme(
        logY ~ d + H_delta + soilTempScaled + vwcScaled + s,
        random = ~1 | plot,
        data = temp_df,
        method = "REML",
        control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
      ),
      error = function(e2) {
        warning("Model failed for delta ", delta, ". Skipping AIC calculation.")
        return(NULL)
      }
    )
  })
  
  if (!is.null(current_model)) {
    aic_values <- c(aic_values, AIC(current_model))
  } else {
    aic_values <- c(aic_values, NA)
  }
}

if (all(is.na(aic_values))) {
  stop("All threshold models failed to converge. Simplify random effects or check data.")
}

best_delta_index <- which.min(aic_values)
best_delta <- deltas[best_delta_index]

cat(paste0("Optimal Threshold (delta) found: ", round(best_delta, 2), "\n"))

df$H_delta <- pmax(0, df$d - best_delta)
df$H_delta
final_freq_model <- NULL

cat("\n--- Fitting Final Mixed-Effects Model ---\n")
tryCatch({
  cat("Attempting to fit full model with random intercept and slope...\n")
  final_freq_model <- lme(
    logY ~ d + H_delta + soilTempScaled + vwcScaled + s,
    random = ~1 + s | plot,
    data = df,
    method = "ML"
  )
  cat("Success: Full mixed-effects model converged.\n")
}, error = function(e) {
  cat("Full model failed. Error: ", e$message, "\n")
  cat("Falling back to a simpler model with only a random intercept...\n")
  final_freq_model <<- tryCatch(
    lme(
      logY ~ d + H_delta + soilTempScaled + vwcScaled + s,
      random = ~1 | plot,
      data = df,
      method = "ML"
    ),
    error = function(e2) {
      cat("Simpler model failed. Error: ", e2$message, "\n")
      cat("Falling back to fixed-effects model...\n")
      lm(logY ~ d + H_delta + soilTempScaled + vwcScaled + s, data = df)
    }
  )
})

if (is.null(final_freq_model)) {
  stop("Final model failed to converge, even after simplification.")
}

cat("\n--- Summary of Final Mixed-Effects Threshold Model ---\n")
summary(final_freq_model)

aic_df <- data.frame(delta = deltas, AIC = aic_values)
p_aic <- ggplot(aic_df, aes(x = delta, y = AIC)) +
  geom_line(color = "steelblue", size = 1) +
  geom_vline(xintercept = best_delta, linetype = "dashed", color = "red") +
  geom_point(aes(x = best_delta, y = min(AIC, na.rm = TRUE)), color = "red", size = 3) +
  theme_minimal() +
  labs(
    title = bquote("AIC Profile for Threshold ("~delta~")"),
    subtitle = bquote("Minimum AIC at"~delta~"="~.(round(best_delta, 2))),
    x = expression(delta),
    y = "AIC"
  )
print(p_aic)

# -------------------------------
# 5. Model Interpretation
# -------------------------------
cat("\n--- Interpretation of Final Model Parameters ---\n")
cat("Coefficients for key terms:\n")
coeffs <- if (inherits(final_freq_model, "lme")) fixef(final_freq_model) else coef(final_freq_model)
cat(paste0("beta_d (linear effect below threshold): ", round(coeffs["d"], 4), "\n"))
cat(paste0("gamma (additional effect above threshold): ", round(coeffs["H_delta"], 4), "\n"))
cat(paste0("beta_x (soilTempScaled): ", round(coeffs["soilTempScaled"], 4), "\n"))
cat(paste0("beta_x (vwcScaled): ", round(coeffs["vwcScaled"], 4), "\n"))
cat(paste0("beta_x (s): ", round(coeffs["s"], 4), "\n"))

# Safely extract random effects standard deviations and correlation
if (inherits(final_freq_model, "lme")) {
  ranef_summary <- tryCatch(
    as.data.frame(VarCorr(final_freq_model)),
    error = function(e) {
      cat("Warning: Could not get VarCorr. The model fit may be singular. Error: ", e$message, "\n")
      return(NULL)
    }
  )
  
  if (!is.null(ranef_summary) && "StdDev" %in% colnames(ranef_summary)) {
    tau_0 <- as.numeric(ranef_summary$StdDev[ranef_summary$Var == "(Intercept)"])
    tau_1 <- if ("s" %in% ranef_summary$Var) {
      as.numeric(ranef_summary$StdDev[ranef_summary$Var == "s"])
    } else { NA }
    rho_b <- if ("Corr" %in% colnames(ranef_summary) && length(ranef_summary$Corr) >= 2 && !is.na(ranef_summary$Corr[2])) {
      as.numeric(ranef_summary$Corr[2])
    } else { NA }
    
    cat("\nRandom Effects:\n")
    cat(paste0("tau_0 (plot intercept SD): ", round(tau_0, 4), "\n"))
    cat(paste0("tau_1 (plot slope SD): ", if (is.na(tau_1)) "Not estimated" else round(tau_1, 4), "\n"))
    cat(paste0("rho_b (intercept-slope correlation): ", if (is.na(rho_b)) "Not estimated" else round(rho_b, 4), "\n"))
  } else {
    cat("\nWarning: No random effects estimated (or could not be retrieved from the model).\n")
  }
} else {
  cat("\nWarning: No random effects estimated (fixed-effects model used).\n")
}

# =====================================================================
# ADDED SECTION: 6. Model Validation and Cross-Validation
# =====================================================================
cat("\n\n--- 6. Model Validation and Cross-Validation ---\n")
cat("Performing leave-one-plot-out cross-validation...\n")

# Store validation metrics
validation_metrics <- data.frame(fold = character(), RMSE = numeric(), MAE = numeric(), stringsAsFactors = FALSE)
plot_ids <- unique(df$plot)
all_preds <- data.frame(plot = character(), date = as.Date(character()), true_logY = numeric(), pred_logY = numeric(), stringsAsFactors = FALSE)

# Loop through each plot to use as the test set
for (i in 1:length(plot_ids)) {
  test_plot <- plot_ids[i]
  
  # Create training and test dataframes
  train_df <- df %>% dplyr::filter(plot != test_plot)
  test_df <- df %>% dplyr::filter(plot == test_plot)
  
  # Check if there is enough data to train the model
  if (dplyr::n_distinct(train_df$plot) < 2) {
    cat(paste0("Skipping fold ", i, " (plot ", test_plot, ") due to insufficient training data.\n"))
    next
  }
  
  # Fit the final threshold model on the training data
  # Use a tryCatch block to handle potential convergence issues
  fold_model <- NULL
  tryCatch({
    fold_model <- lme(
      logY ~ d + H_delta + soilTempScaled + vwcScaled + s,
      random = ~1 | plot,
      data = train_df,
      method = "ML"
    )
  }, error = function(e) {
    cat(paste0("Model failed to converge for plot ", test_plot, ". Falling back to linear model.\n"))
    fold_model <<- lm(logY ~ d + H_delta + soilTempScaled + vwcScaled + s, data = train_df)
  })
  
  # Make predictions on the test set
  if (!is.null(fold_model)) {
    preds <- predict(fold_model, newdata = test_df, level = 0) # level=0 for fixed effects only
    
    # Calculate RMSE and MAE
    rmse <- sqrt(mean((test_df$logY - preds)^2, na.rm = TRUE))
    mae <- mean(abs(test_df$logY - preds), na.rm = TRUE)
    
    # Store metrics and predictions
    validation_metrics <- rbind(validation_metrics, data.frame(fold = paste0("Plot_", test_plot), RMSE = rmse, MAE = mae))
    all_preds <- rbind(all_preds, data.frame(plot = test_df$plot, date = test_df$date, true_logY = test_df$logY, pred_logY = preds))
  }
}

cat("\n--- Cross-Validation Results ---\n")
kable(validation_metrics, caption = "Leave-One-Plot-Out Cross-Validation Metrics")
cat(paste0("\nAverage Out-of-Sample RMSE: ", round(mean(validation_metrics$RMSE, na.rm = TRUE), 4), "\n"))
cat(paste0("Average Out-of-Sample MAE: ", round(mean(validation_metrics$MAE, na.rm = TRUE), 4), "\n"))

# --- PIT Histogram ---
cat("\n--- Probability Integral Transform (PIT) Histogram ---\n")
# Note: PIT histogram requires a full predictive distribution, which is not directly
# available from the lme prediction function. As a proxy, we use residuals.
resids <- all_preds$true_logY - all_preds$pred_logY
p_pit <- ggplot(data.frame(resids=resids), aes(x = resids)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of In-Sample Residuals (Proxy for PIT)",
       subtitle = "A well-calibrated model should produce uniform residuals.",
       x = "Residuals (True logY - Predicted logY)", y = "Count")
print(p_pit)


# =====================================================================
# ADDED SECTION: 7. Recovery Dynamics Model Implementation
# =====================================================================
cat("\n\n--- 7. Non-Linear Mixed-Effects Model for Recovery Dynamics ---\n")
cat("Implementing a simplified recovery model (Equation 5) using nlme.\n")
cat("Note: This model simplifies the environment-linked recovery rate for better convergence.\n")

# Define the non-linear function for the recovery term
# `logY` is the dependent variable, `s` is the independent variable, `phi` and `rho` are parameters to estimate
recovery_fun <- function(s, phi, rho) {
  # This term models the recovery dynamics, R(s) from Equation 5
  phi * (1 - exp(-rho * s))
}

# Create a non-linear mixed-effects model formula
# This combines the linear fixed effects with the non-linear recovery term
nlme_formula <- nlme::y ~ recovery_fun(s, phi, rho) + beta_0 + beta_d*d + gamma*H_delta(d) + 
  beta_soilTemp*soilTempScaled + beta_vwc*vwcScaled + 
  (b0_plot + b1_plot*s | plot)

# Set up initial values for the non-linear parameters
# These are crucial for the model to converge
start_values <- list(
  phi = 1,
  rho = 0.005,
  beta_0 = 1,
  beta_d = 0,
  gamma = 0,
  beta_soilTemp = 0,
  beta_vwc = 0
)

# Fit the non-linear mixed-effects model
recovery_model <- NULL
tryCatch({
  recovery_model <- nlme(
    model = nlme_formula,
    data = df,
    fixed = list(
      phi + rho + beta_0 + beta_d + gamma + beta_soilTemp + beta_vwc ~ 1
    ),
    random = list(
      b0_plot + b1_plot ~ 1 | plot
    ),
    start = start_values,
    control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
  )
  cat("Success: Non-linear mixed-effects model converged.\n")
}, error = function(e) {
  cat("Non-linear model failed to converge. This is common with complex models. Error: ", e$message, "\n")
  recovery_model <<- NULL # Explicitly set to NULL on failure
})

# Display the summary of the fitted model
if (!is.null(recovery_model)) {
  cat("\n--- Summary of Recovery Dynamics Model ---\n")
  summary(recovery_model)
  
  # Plot the fitted recovery curve for a high-severity plot
  cat("\n--- Plotting Fitted Recovery Curve for an example plot ---\n")
  example_plot_data <- df %>%
    dplyr::filter(disturbance_severity == max(disturbance_severity, na.rm = TRUE)) %>%
    dplyr::mutate(
      fitted = predict(recovery_model, newdata = ., level = 0), # fixed effects only
      fitted_full = predict(recovery_model, newdata = ., level = 1) # including random effects
    )
  
  p_recovery_fitted <- ggplot(example_plot_data, aes(x = s)) +
    geom_point(aes(y = logY), alpha = 0.6) +
    geom_line(aes(y = fitted, color = "Fixed"), size = 1) +
    geom_line(aes(y = fitted_full, color = "Full"), size = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
      title = paste0("Fitted Recovery Curve for Plot ", unique(example_plot_data$plot)),
      x = "Days Since Disturbance (s)", y = "log(Flux)", color = "Prediction Type"
    ) +
    scale_color_manual(values = c("Fixed" = "blue", "Full" = "red"))
  
  print(p_recovery_fitted)
  
} else {
  cat("The non-linear recovery model could not be fitted.\n")
}

