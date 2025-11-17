### Occupancy Model Code ###

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, tidyr, readr, nasapower, lubridate, unmarked, terra, sf, MuMIn, tidyverse,
  ggplot2, hms, purrr, corrr, patchwork, cowplot, ordinal, emmeans, AICcmodavg,
  readxl, ggthemes, rnaturalearth, ggspatial
)

# Load cleaned data
survey_data_wide <- read.csv("data/occupancy_model_data_cleaned.csv")

sig_cols <- survey_data_wide %>% select(sig_abundance1:sig_abundance6)
get_max_ordinal <- function(x) {
  valid <- x[!is.na(x) & x %in% 0:4]
  if (length(valid) == 0) return(NA) else return(max(valid))
}
survey_data_wide$sig_abundance_max_cat <- apply(sig_cols, 1, get_max_ordinal)

z.transform <- function(x) {  # defined early because we use it here
  (x - mean(x, na.rm = TRUE)) / (2 * sd(x, na.rm = TRUE))
}
survey_data_wide$sig_abundance_max_cat_z <- z.transform(survey_data_wide$sig_abundance_max_cat)  # <<< AICc patch keeps

### Detection matrix
y_matrix <- as.matrix(survey_data_wide %>%
                        select(bib_present1, bib_present2, bib_present3,
                               bib_present4, bib_present5, bib_present6))

# Site-level covariates (standardised)
survey_data_wide$elevation_std <- z.transform(survey_data_wide$elevation)
survey_data_wide$canopy_cover_site_std <- z.transform(survey_data_wide$canopy_cover_site)
survey_data_wide$fire_severity_std <- z.transform(survey_data_wide$fire_severity)

# Observation-level covariates (standardise)
for (i in 1:6) {
  survey_data_wide[[paste0("daily_temp", i, "_std")]]      <- z.transform(survey_data_wide[[paste0("daily_temp", i)]])
  survey_data_wide[[paste0("daily_precip", i, "_std")]]    <- z.transform(survey_data_wide[[paste0("daily_precip", i)]])
  survey_data_wide[[paste0("weekly_precip", i, "_std")]]   <- z.transform(survey_data_wide[[paste0("weekly_precip", i)]])
  survey_data_wide[[paste0("next_day_precip", i, "_std")]] <- z.transform(survey_data_wide[[paste0("next_day_precip", i)]])
  survey_data_wide[[paste0("wind", i, "_std")]]            <- z.transform(survey_data_wide[[paste0("wind", i)]])
  survey_data_wide[[paste0("humidity", i, "_std")]]        <- z.transform(survey_data_wide[[paste0("humidity", i)]])
  survey_data_wide[[paste0("moisture_level", i, "_std")]]  <- z.transform(survey_data_wide[[paste0("moisture_level", i)]])
  survey_data_wide[[paste0("cloud_cover", i, "_std")]]     <- z.transform(survey_data_wide[[paste0("cloud_cover", i)]])
}

# ----- use numeric z variable in site_covs -----
site_covs <- survey_data_wide %>% select(elevation_std, sig_abundance_max_cat_z)

### DAILY PRECIP model
obs_covs_daily <- list(
  daily_temp = as.matrix(survey_data_wide[, paste0("daily_temp", 1:6, "_std")]),
  daily_precip = as.matrix(survey_data_wide[, paste0("daily_precip", 1:6, "_std")]),
  wind = as.matrix(survey_data_wide[, paste0("wind", 1:6, "_std")]),
  humidity = as.matrix(survey_data_wide[, paste0("humidity", 1:6, "_std")]),
  moisture_level = as.matrix(survey_data_wide[, paste0("moisture_level", 1:6, "_std")]),
  cloud_cover = as.matrix(survey_data_wide[, paste0("cloud_cover", 1:6, "_std")]),
  cosHour = as.matrix(survey_data_wide[, paste0("cosHour", 1:6)]),
  sinHour = as.matrix(survey_data_wide[, paste0("sinHour", 1:6)])
)

umf_daily <- unmarkedFrameOccu(
  y = y_matrix,
  siteCovs = site_covs,
  obsCovs = obs_covs_daily
)


fit_daily <- occu(~ daily_temp + daily_precip + wind + humidity + moisture_level + cloud_cover + cosHour + sinHour ~
                    elevation_std + sig_abundance_max_cat_z, data = umf_daily)

### WEEKLY PRECIP model
obs_covs_weekly <- list(
  daily_temp = as.matrix(survey_data_wide[, paste0("daily_temp", 1:6, "_std")]),
  weekly_precip = as.matrix(survey_data_wide[, paste0("weekly_precip", 1:6, "_std")]),
  wind = as.matrix(survey_data_wide[, paste0("wind", 1:6, "_std")]),
  humidity = as.matrix(survey_data_wide[, paste0("humidity", 1:6, "_std")]),
  moisture_level = as.matrix(survey_data_wide[, paste0("moisture_level", 1:6, "_std")]),
  cloud_cover = as.matrix(survey_data_wide[, paste0("cloud_cover", 1:6, "_std")]),
  cosHour = as.matrix(survey_data_wide[, paste0("cosHour", 1:6)]),
  sinHour = as.matrix(survey_data_wide[, paste0("sinHour", 1:6)])
)

umf_weekly <- unmarkedFrameOccu(
  y = y_matrix,
  siteCovs = site_covs,
  obsCovs = obs_covs_weekly
)


fit_weekly <- occu(~ daily_temp + weekly_precip + wind + humidity + moisture_level + cloud_cover + cosHour + sinHour ~
                     elevation_std + sig_abundance_max_cat_z, data = umf_weekly)

### NEXT DAY PRECIP model
obs_covs_nextday <- list(
  daily_temp = as.matrix(survey_data_wide[, paste0("daily_temp", 1:6, "_std")]),
  next_day_precip = as.matrix(survey_data_wide[, paste0("next_day_precip", 1:6, "_std")]),
  wind = as.matrix(survey_data_wide[, paste0("wind", 1:6, "_std")]),
  humidity = as.matrix(survey_data_wide[, paste0("humidity", 1:6, "_std")]),
  moisture_level = as.matrix(survey_data_wide[, paste0("moisture_level", 1:6, "_std")]),
  cloud_cover = as.matrix(survey_data_wide[, paste0("cloud_cover", 1:6, "_std")]),
  cosHour = as.matrix(survey_data_wide[, paste0("cosHour", 1:6)]),
  sinHour = as.matrix(survey_data_wide[, paste0("sinHour", 1:6)])
)

umf_nextday <- unmarkedFrameOccu(
  y = y_matrix,
  siteCovs = site_covs,
  obsCovs = obs_covs_nextday
)


fit_nextday <- occu(~ daily_temp + next_day_precip + wind + humidity + moisture_level + cloud_cover + cosHour + sinHour ~
                      elevation_std + sig_abundance_max_cat_z, data = umf_nextday)

### Compare AICc 
aic_vals <- c(
  daily    = MuMIn::AICc(fit_daily),
  weekly   = MuMIn::AICc(fit_weekly),
  next_day = MuMIn::AICc(fit_nextday)
)
print(round(aic_vals, 2))

### USE NEXT DAY AS LOWEST AICc


##############################
# Occupancy model #
##############################

library(dplyr)
library(unmarked)
library(MuMIn)

options(na.action = "na.fail") 

#---------------------------
# Load & prepare data
#---------------------------
survey_data_wide <- read.csv("data/occupancy_model_data_cleaned.csv")

# Max Crinia category across visits (0–4) 
sig_cols <- survey_data_wide %>% dplyr::select(sig_abundance1:sig_abundance6)
get_max_ordinal <- function(x) {
  valid <- x[!is.na(x) & x %in% 0:4]
  if (length(valid) == 0) return(NA_integer_) else return(max(valid))
}
survey_data_wide$sig_abundance_max_cat <- apply(sig_cols, 1, get_max_ordinal)


z.transform <- function(x) (x - mean(x, na.rm = TRUE)) / (2 * stats::sd(x, na.rm = TRUE))
survey_data_wide$sig_abundance_max_cat_z <- z.transform(survey_data_wide$sig_abundance_max_cat)

# Response matrix
y <- survey_data_wide %>%
  dplyr::select(dplyr::starts_with("bib_present")) %>%
  as.matrix()

# Site-level covariates (std where appropriate)
survey_data_wide$elevation_std         <- z.transform(survey_data_wide$elevation)
survey_data_wide$canopy_cover_site_std <- z.transform(survey_data_wide$canopy_cover_site)
survey_data_wide$fire_severity_std     <- z.transform(survey_data_wide$fire_severity)


site_covs <- survey_data_wide %>%
  dplyr::select(elevation_std, canopy_cover_site_std, fire_severity_std, sig_abundance_max_cat_z)

# Observation-level covariates (std where appropriate)
for (i in 1:6) {
  survey_data_wide[[paste0("daily_temp", i, "_std")]]      <- z.transform(survey_data_wide[[paste0("daily_temp", i)]])
  survey_data_wide[[paste0("next_day_precip", i, "_std")]] <- z.transform(survey_data_wide[[paste0("next_day_precip", i)]])
  survey_data_wide[[paste0("wind", i, "_std")]]            <- z.transform(survey_data_wide[[paste0("wind", i)]])
  survey_data_wide[[paste0("humidity", i, "_std")]]        <- z.transform(survey_data_wide[[paste0("humidity", i)]])
  survey_data_wide[[paste0("moisture_level", i, "_std")]]  <- z.transform(survey_data_wide[[paste0("moisture_level", i)]])
  survey_data_wide[[paste0("cloud_cover", i, "_std")]]     <- z.transform(survey_data_wide[[paste0("cloud_cover", i)]])
}

obs_covs <- list(
  daily_temp          = as.matrix(survey_data_wide[, paste0("daily_temp", 1:6, "_std")]),
  next_day_precip     = as.matrix(survey_data_wide[, paste0("next_day_precip", 1:6, "_std")]),
  wind                = as.matrix(survey_data_wide[, paste0("wind", 1:6, "_std")]),
  humidity            = as.matrix(survey_data_wide[, paste0("humidity", 1:6, "_std")]),
  moisture_level      = as.matrix(survey_data_wide[, paste0("moisture_level", 1:6, "_std")]),
  cloud_cover         = as.matrix(survey_data_wide[, paste0("cloud_cover", 1:6, "_std")]),
  cosHour             = as.matrix(survey_data_wide[, paste0("cosHour", 1:6)]),
  sinHour             = as.matrix(survey_data_wide[, paste0("sinHour", 1:6)]),
  habitat_disturbance = as.matrix(survey_data_wide[, paste0("habitat_disturbance", 1:6)])
)

# Unmarked frame
umf <- unmarkedFrameOccu(y = y, siteCovs = site_covs, obsCovs = obs_covs)

#---------------------------
# Candidate set
#---------------------------
# Detection formula: same structure across all models
det_formula <- paste("~", paste(names(obs_covs), collapse = " + "))


occ_covs <- c("elevation_std", "sig_abundance_max_cat_z", "fire_severity_std", "canopy_cover_site_std")
occ_formulas <- unlist(lapply(1:length(occ_covs), function(i) {
  combn(occ_covs, i, FUN = function(x) paste(x, collapse = " + "), simplify = TRUE)
}))

model_list <- lapply(occ_formulas, function(occ) {
  f <- as.formula(paste(det_formula, "~", occ))
  tryCatch(occu(f, data = umf), error = function(e) NULL)
})
model_list <- model_list[!vapply(model_list, is.null, logical(1))]
names(model_list) <- occ_formulas[!vapply(model_list, is.null, logical(1))]

#---------------------------
# Model selection (AICc)
#---------------------------
fit_list <- fitList(fits = model_list)

ms_unmarked <- MuMIn::model.sel(model_list, rank = "AICc")
print(ms_unmarked)



#---------------------------
# Model averaging (ΔAICc < 2)
#---------------------------
ms_tab <- MuMIn::model.sel(model_list, rank = "AICc")


top_models <- MuMIn::get.models(ms_tab, subset = delta < 2)  
summary(top_models)

top_models <- top_models[1:2]
avg_results <- MuMIn::model.avg(top_models, rank = "AICc")
summary(avg_results)
confint(avg_results)

# ---- Collinearity checks: elevation vs fire severity ----
sc <- siteCovs(umf)

cor_pearson  <- cor(sc$elevation_std, sc$fire_severity_std, use = "complete.obs", method = "pearson")
cor_spearman <- cor(sc$elevation_std, sc$fire_severity_std, use = "complete.obs", method = "spearman")
cat("Pearson r (elev vs fire):", round(cor_pearson, 3), "\n")
cat("Spearman rho (elev vs fire):", round(cor_spearman, 3), "\n")

ct <- cor.test(sc$elevation_std, sc$fire_severity_std, method = "pearson")
print(ct)

library(ggplot2)
ggplot(sc, aes(x = elevation_std, y = fire_severity_std)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Elevation (standardised)", y = "Fire severity (standardised)")

# ----- VIF -----
set.seed(1)
dummy <- rnorm(nrow(sc))
m_vif <- lm(dummy ~ elevation_std + fire_severity_std + canopy_cover_site_std + sig_abundance_max_cat_z, data = sc)
vifs <- car::vif(m_vif)
print(vifs)
# Rule-of-thumb: |r| > 0.7 high; VIF > 5 (or >10) concerning.

####################################
# VISUALISE THE PREDICTIONS
####################################
library(unmarked)
library(ggplot2)
library(patchwork)

## back-transform standardized variables (works with scale())
unscale_vec <- function(z, scaled_vec) {
  cen <- attr(scaled_vec, "scaled:center")
  sca <- attr(scaled_vec, "scaled:scale")
  if (!is.null(cen) && !is.null(sca)) z * sca + cen else z
}

## =================================================================================
## 1) FIT THE AVERAGED OCCUPANCY MODEL
##    ψ ~ elevation_std + fire_severity_std + sig_abundance_max_cat_z
##    p ~ daily_temp + next_day_precip + wind + humidity +
##         moisture_level + cloud_cover + habitat_disturbance + cosHour + sinHour
## =================================================================================
refit_model <- occu(
  ~ daily_temp + next_day_precip + wind + humidity +
    moisture_level + cloud_cover + habitat_disturbance + cosHour + sinHour ~
    elevation_std + fire_severity_std + sig_abundance_max_cat_z,
  data = umf
)

## =================================================================================
## 2) MARGINAL EFFECT: OCCUPANCY vs ELEVATION (fire=0, Crinia_z=0)
## =================================================================================
elev_center <- mean(survey_data_wide$elevation, na.rm = TRUE)
elev_scale  <- 2 * sd(survey_data_wide$elevation, na.rm = TRUE)

elev_raw_seq <- seq(
  from = min(survey_data_wide$elevation, na.rm = TRUE),
  to   = max(survey_data_wide$elevation, na.rm = TRUE),
  length.out = 200
)
elev_seq <- (elev_raw_seq - elev_center) / elev_scale
stopifnot(length(elev_seq) == 200, all(is.finite(elev_seq)))

pred_data_elev <- data.frame(
  elevation_std            = elev_seq,
  fire_severity_std        = 0,
  sig_abundance_max_cat_z  = 0
)
pred_elev <- predict(refit_model, type = "state", newdata = pred_data_elev)

df_elev <- data.frame(
  elevation_raw = elev_raw_seq,
  estimate = pred_elev$Predicted,
  lower    = pred_elev$lower,
  upper    = pred_elev$upper
)

p_occ_elev <- ggplot(df_elev, aes(x = elevation_raw, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.4) +
  geom_line(linewidth = 1) +
  scale_x_continuous(limits = c(0, 1700), breaks = seq(0, 1700, by = 500)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Elevation (m)", y = "Occupancy probability") +
  theme_classic()

## =================================================================================
## 3) MARGINAL EFFECT: OCCUPANCY vs FIRE SEVERITY (elev=0, Crinia_z=0)
## =================================================================================
fire_seq <- seq(min(siteCovs(umf)$fire_severity_std, na.rm = TRUE),
                max(siteCovs(umf)$fire_severity_std, na.rm = TRUE),
                length.out = 200)

pred_data_fire <- data.frame(
  elevation_std            = 0,
  fire_severity_std        = fire_seq,
  sig_abundance_max_cat_z  = 0
)
pred_fire <- predict(refit_model, type = "state", newdata = pred_data_fire)

fire_center <- mean(survey_data_wide$fire_severity, na.rm = TRUE)
fire_scale  <- 2 * sd(survey_data_wide$fire_severity, na.rm = TRUE)
fire_raw    <- fire_seq * fire_scale + fire_center

df_fire <- data.frame(
  fire_raw = fire_raw,
  estimate = pred_fire$Predicted,
  lower = pred_fire$lower,
  upper = pred_fire$upper
)

p_occ_fire <- ggplot(df_fire, aes(x = fire_raw, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "tomato", alpha = 0.2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Fire severity (GEEBAM index)",
    y ="Occupancy probability "
  ) +
  theme_classic()

## =================================================================================
## 4) MARGINAL EFFECT: OCCUPANCY vs CRINIA (sig_abundance_max_cat_z)
##    Hold elevation at mean (0), fire at mean (0)
## =================================================================================
crinia_z_seq <- seq(
  from = min(siteCovs(umf)$sig_abundance_max_cat_z, na.rm = TRUE),
  to   = max(siteCovs(umf)$sig_abundance_max_cat_z, na.rm = TRUE),
  length.out = 200
)

pred_data_crinia <- data.frame(
  elevation_std            = 0,
  fire_severity_std        = 0,
  sig_abundance_max_cat_z  = crinia_z_seq
)

pred_crinia <- predict(refit_model, type = "state", newdata = pred_data_crinia)

# Back-transform z to the original 0–4 scale
crinia_center <- mean(survey_data_wide$sig_abundance_max_cat, na.rm = TRUE)
crinia_scale  <- 2 * sd(survey_data_wide$sig_abundance_max_cat, na.rm = TRUE)
crinia_raw    <- crinia_z_seq * crinia_scale + crinia_center

df_crinia <- data.frame(
  crinia_raw = crinia_raw,
  estimate   = pred_crinia$Predicted,
  lower      = pred_crinia$lower,
  upper      = pred_crinia$upper
)

## --- anchor points exactly at the five categories (0,1,2,3,4) ---
cat_levels <- 0:4
crinia_z_at_levels <- (cat_levels - crinia_center) / crinia_scale
pred_levels <- predict(
  refit_model, type = "state",
  newdata = data.frame(
    elevation_std = 0,
    fire_severity_std = 0,
    sig_abundance_max_cat_z = crinia_z_at_levels
  )
)

pts_crinia <- data.frame(
  crinia_raw = cat_levels,
  estimate   = pred_levels$Predicted,
  lower      = pred_levels$lower,
  upper      = pred_levels$upper
)

## plot with categorical tick labels
p_occ_crinia <- ggplot(df_crinia, aes(x = crinia_raw, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "darkseagreen3", alpha = 0.25) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = 0:4,
    labels = c("0", "1–5", "6–20", "21–50", "51–100"),
    limits = c(0, 4)
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = expression(italic(Crinia~signifera) ~ "maximum count category (individuals)"),
    y = "Occupancy probability"
  ) + 
  theme_classic()

p_occ_crinia
## =================================================================================
## 5) DETECTION vs TIME-OF-DAY (cosHour / sinHour)
## =================================================================================
pred_df <- data.frame(hour = seq(0, 23.99, length.out = 400))
pred_df$cosHour <- cos((2 * pi * pred_df$hour) / 24)
pred_df$sinHour <- sin((2 * pi * pred_df$hour) / 24)

det_baseline <- data.frame(
  daily_temp          = 0,
  next_day_precip     = 0,
  wind                = 0,
  humidity            = 0,
  moisture_level      = 0,
  cloud_cover         = 0,
  habitat_disturbance = 0
)
det_baseline <- det_baseline[rep(1, nrow(pred_df)), ]
pred_det_df <- cbind(pred_df, det_baseline)

pred_det <- predict(refit_model, type = "det", newdata = pred_det_df)

df_det <- data.frame(
  hour  = pred_df$hour,
  pred  = pred_det$Predicted,
  lower = pred_det$lower,
  upper = pred_det$upper
)

p_det <- ggplot(df_det, aes(x = hour, y = pred)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "coral3", alpha = 0.3) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 24, by = 4),
                     labels = c("12AM","4AM","8AM","12PM","4PM","8PM","12AM")) +
  labs(
    x = "Time of day",
    y = "Detection probability") +
  theme_classic()

## =================================================================================
## SAVE
## =================================================================================
library(AICcmodavg)
gof <- mb.gof.test(refit_model, nsim = 1000)
gof$chat  # if ~1.5–2+, expect wider SEs
print(gof)
### they fit fine P value - 0.394 and no over dispersion (c-hat 0.88)

# Marginal effects (single Crinia level → now mean z=0)
p_combined <- p_occ_elev / p_occ_crinia + plot_annotation(tag_levels = "A")
p_combined
ggsave("bibronii_occupancy_vs_crinia.png", p_combined,
       width = 7, height = 8, dpi = 300)

# Detection vs time
p_det
ggsave("bibronii_detection_vs_time.png", p_det, width = 8, height = 5, dpi = 300)

library(ggplot2)

# =========================
# Crinia × Elevation figure (with point legend)
# =========================
library(dplyr)
library(ggplot2)

# --- Centers/scales used for Gelman transform ---
elev_center  <- mean(survey_data_wide$elevation, na.rm = TRUE)
elev_scale   <- 2 * sd(survey_data_wide$elevation, na.rm = TRUE)

crinia_center <- mean(survey_data_wide$sig_abundance_max_cat, na.rm = TRUE)
crinia_scale  <- 2 * sd(survey_data_wide$sig_abundance_max_cat, na.rm = TRUE)

elev_raw_vals <- as.numeric(quantile(survey_data_wide$elevation,
                                     probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
elev_z_vals   <- (elev_raw_vals - elev_center) / elev_scale

crinia_z_seq <- seq(
  from = min(siteCovs(umf)$sig_abundance_max_cat_z, na.rm = TRUE),
  to   = max(siteCovs(umf)$sig_abundance_max_cat_z, na.rm = TRUE),
  length.out = 200
)

# Build prediction grid (fire at mean = 0)
grid <- expand.grid(
  elevation_std            = elev_z_vals,
  sig_abundance_max_cat_z  = crinia_z_seq,
  fire_severity_std        = 0
)

pred <- predict(refit_model, type = "state", newdata = grid)

df_smooth <- cbind(grid, pred) |>
  dplyr::transmute(
    elevation_std,
    crinia_raw = sig_abundance_max_cat_z * crinia_scale + crinia_center,  # back-transform to 0–4
    Predicted, lower, upper
  )

# Label the three elevation groups for the lines
df_smooth$elev_lab <- factor(
  elev_raw_vals[match(df_smooth$elevation_std, elev_z_vals)],
  levels = elev_raw_vals,
  labels = c("Low elevation", "Mid elevation", "High elevation")
)

# Elevation quartile groups collapsed to Low/Mid/High
cuts <- quantile(survey_data_wide$elevation, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
elev_grp_raw <- cut(survey_data_wide$elevation, breaks = cuts, include.lowest = TRUE,
                    labels = c("Low (Q1)","Q1–Q2","Q2–Q3","High (Q3)"))

site_pred <- predict(refit_model, type = "state", newdata = siteCovs(umf))

site_points <- data.frame(
  psi        = site_pred$Predicted,
  crinia_x   = as.numeric(survey_data_wide$sig_abundance_max_cat),  # 0..4 on the x-axis
  elev_lab   = dplyr::case_when(
    elev_grp_raw == "Low (Q1)"                 ~ "Low elevation",
    elev_grp_raw %in% c("Q1–Q2","Q2–Q3")       ~ "Mid elevation",
    elev_grp_raw == "High (Q3)"                ~ "High elevation",
    TRUE                                       ~ NA_character_
  )
) |>
  dplyr::filter(!is.na(elev_lab))

site_points$elev_lab <- factor(site_points$elev_lab,
                               levels = c("Low elevation","Mid elevation","High elevation"))

# Plot (with separate Type legend for points vs line)
p_crinia_by_elev <- ggplot(df_smooth, aes(x = crinia_raw, y = Predicted, group = elev_lab)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = elev_lab),
              alpha = 0.20, colour = NA, show.legend = FALSE) +
  geom_line(aes(colour = elev_lab, linetype = "Model curve"), linewidth = 1) +  # adds "Type" legend entry
  # points = site-level predicted ψ at each site's observed Crinia category
  geom_jitter(
    data = site_points,
    aes(x = crinia_x, y = psi, colour = elev_lab, shape = "Sites (predicted \u03C8)"),
    width = 0.08, alpha = 0.35, size = 1.8, inherit.aes = FALSE
  ) +
  scale_x_continuous(
    breaks = 0:4,
    labels = c("0","1–5","6–20","21–50","51–100"),
    limits = c(0, 4)
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_linetype_manual(name = "Type", values = 1) +     # "Model curve"
  scale_shape_manual(name = "Type", values = 16) +       # "Sites (predicted ψ)"
  guides(
    colour  = guide_legend(title = "Elevation", order = 1),
    linetype= guide_legend(order = 2),
    shape   = guide_legend(order = 2)
  ) +
  labs(
    x = expression(italic(Crinia~signifera)~" maximum count category (individuals)"),
    y = "Predicted occupancy probability"
  ) +
  theme_classic()

p_crinia_by_elev

ggsave("bibronii_crinia_by_elevation_with_point_legend.png", p_crinia_by_elev, width = 7, height = 5, dpi = 300)
