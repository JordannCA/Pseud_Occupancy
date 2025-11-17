
###### CLMM MODEL FOR BIBRONII ABUNDANCE ######
library(dplyr)
library(tidyr)
library(ordinal)
library(readr)
library(AICcmodavg)
library(broom.mixed)
library(lubridate)
library(MuMIn)   

options(na.action = "na.fail")  # MuMIn best practice

z.transform <- function(x) (x - mean(x, na.rm = TRUE)) / (2 * stats::sd(x, na.rm = TRUE))

df <- read_csv("data/occupancy_model_data_cleaned.csv")
names(df)

# Pivot to long format for survey-level variables
df_long <- df %>%
  pivot_longer(
    cols = c(starts_with("bib_abundance"), starts_with("sig_abundance"),
             starts_with("daily_temp"), starts_with("time"), starts_with("date"),
             starts_with("next_day_precip")),
    names_to = c(".value", "survey_number"),
    names_pattern = "(.*?)(\\d+)"
  ) %>%
  filter(!is.na(bib_abundance)) %>%
  mutate(
    survey_number       = as.integer(survey_number),
    bib_abundance_clean = pmin(round(bib_abundance), 4),
    bib_abundance_cat   = factor(
      bib_abundance_clean,
      levels = 0:4,
      labels = c("0", "1–5", "6–20", "21–50", "51–100"),
      ordered = TRUE
    )
  )

# Day vs night
library(lubridate)
morning <- hm("07:00")
evening <- hm("19:00")
df_long$night <- ifelse(df_long$time > morning & df_long$time < evening, 0, 1)

# Keep surveys with a defined count category for bibs
cats <- c(0,1,2,3,4)
df_long <- df_long[df_long$bib_abundance %in% cats, ]
nrow(df_long)

# ---- site-level max Crinia category (0–4) ----
max.sig.func <- function(site){
  v <- df_long$sig_abundance[df_long$site_ID == site]
  v <- v[!is.na(v) & v %in% 0:4]
  if (length(v) == 0) return(NA_integer_) else return(max(v))
}
df_long$max_sig_abundance <- NA_integer_
for(i in 1:nrow(df_long)){
  df_long$max_sig_abundance[i] <- max.sig.func(site = df_long$site_ID[i])
}

z.transform <- function(x) {
  x <- as.numeric(x)
  m <- mean(x, na.rm = TRUE); s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / (2 * s)
}

# Apply to all continuous covariates
df_long <- df_long %>%
  mutate(
    elevation_std        = z.transform(elevation),
    canopy_std           = z.transform(canopy_cover_site),
    fire_severity_std    = z.transform(fire_severity),
    sig_abundance_std    = z.transform(max_sig_abundance),  # site-level 0–4
    temp_std             = z.transform(daily_temp),
    next_day_precip_std  = z.transform(next_day_precip)
  )

table(df_long$bib_abundance)

# Presence-only for abundance model
df_abund <- df_long %>%
  filter(bib_abundance_cat != "0") %>%
  mutate(
    site_ID = as.factor(site_ID),
    bib_abundance_cat = droplevels(bib_abundance_cat)
  )

table(df_abund$bib_abundance_cat)

# -------------------------------
# Candidate CLMMs
# -------------------------------
fit_clmm <- function(formula) clmm(formula, data = df_abund)

m1  <- fit_clmm(bib_abundance_cat ~ elevation_std + (1 | site_ID))                                   # Elevation only
m2  <- fit_clmm(bib_abundance_cat ~ elevation_std + fire_severity_std + (1 | site_ID))               # + Fire
m3  <- fit_clmm(bib_abundance_cat ~ elevation_std * sig_abundance_std + (1 | site_ID))               # Elev × Crinia
m4  <- fit_clmm(bib_abundance_cat ~ elevation_std * sig_abundance_std + fire_severity_std + (1 | site_ID)) # + Fire
m5  <- fit_clmm(bib_abundance_cat ~ elevation_std + night + (1 | site_ID))                           # + Night
m6  <- fit_clmm(bib_abundance_cat ~ elevation_std + night + fire_severity_std + (1 | site_ID))       # + Night + Fire
m7  <- fit_clmm(bib_abundance_cat ~ elevation_std + next_day_precip_std + (1 | site_ID))             # + Rain
m8  <- fit_clmm(bib_abundance_cat ~ elevation_std + next_day_precip_std + fire_severity_std + (1 | site_ID)) # + Rain + Fire
m9  <- fit_clmm(bib_abundance_cat ~ elevation_std * sig_abundance_std + night + (1 | site_ID))       # Elev×Crinia + Night
m10 <- fit_clmm(bib_abundance_cat ~ elevation_std * sig_abundance_std + night + fire_severity_std + (1 | site_ID)) # + Fire
m11 <- fit_clmm(bib_abundance_cat ~ elevation_std * sig_abundance_std + next_day_precip_std + (1 | site_ID))       # Elev×Crinia + Rain
m12 <- fit_clmm(bib_abundance_cat ~ elevation_std * sig_abundance_std + next_day_precip_std + fire_severity_std + (1 | site_ID)) # + Fire
m13 <- fit_clmm(bib_abundance_cat ~ elevation_std + night + next_day_precip_std + (1 | site_ID))     # Elev + Night + Rain
m14 <- fit_clmm(bib_abundance_cat ~ elevation_std + night + next_day_precip_std + fire_severity_std + (1 | site_ID)) # + Fire
m15 <- fit_clmm(bib_abundance_cat ~ elevation_std * sig_abundance_std + night + next_day_precip_std + (1 | site_ID))  # Elev×Crinia + Night + Rain
m16 <- fit_clmm(bib_abundance_cat ~ elevation_std * sig_abundance_std + night + next_day_precip_std + fire_severity_std + (1 | site_ID)) # + Fire

model_list  <- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)
model_names <- c(
  "Elevation only",
  "Elevation + Fire severity",
  "Elevation × Crinia",
  "Elevation × Crinia + Fire severity",
  "Elevation + Night",
  "Elevation + Night + Fire severity",
  "Elevation + Next-day rain",
  "Elevation + Next-day rain + Fire severity",
  "Elevation × Crinia + Night",
  "Elevation × Crinia + Night + Fire severity",
  "Elevation × Crinia + Next-day rain",
  "Elevation × Crinia + Next-day rain + Fire severity",
  "Elevation + Night + Next-day rain",
  "Elevation + Night + Next-day rain + Fire severity",
  "Elevation × Crinia + Night + Next-day rain",
  "Elevation × Crinia + Night + Next-day rain + Fire severity"
)
names(model_list) <- model_names   

# -------------------------------
# AICc table 
# -------------------------------
aic_tab <- AICcmodavg::aictab(
  cand.set   = model_list,
  modnames   = model_names,
  second.ord = TRUE     # <<< AICc (small-sample corrected)
)
print(aic_tab)

# -------------------------------
# Model selection & averaging (AICc)
# -------------------------------
ms <- MuMIn::model.sel(model_list, rank = "AICc")     
print(ms)

top_models <- MuMIn::get.models(ms, subset = delta < 2)
avg_results <- MuMIn::model.avg(top_models, rank = "AICc")  
summary(avg_results)

# tables 
Table_S6 <- data.frame(
  Model  = rownames(ms),
  K      = ms$df,
  AICc   = round(ms$AICc, 3),   
  Delta  = round(ms$delta, 3),
  Weight = round(ms$weight, 3)
)
Table_S6 <- Table_S6[order(Table_S6$AICc), ]
print(Table_S6, row.names = FALSE)

# Table S7
sum_avg  <- summary(avg_results)
coef_cond <- as.data.frame(sum_avg$coefmat.subset)
coef_cond$Parameter <- rownames(coef_cond)
coef_cond <- coef_cond[!grepl("\\|", coef_cond$Parameter), ]
Table_S7 <- coef_cond |>
  dplyr::transmute(
    Parameter,
    Estimate = round(Estimate, 3),
    SE       = round(`Std. Error`, 3),
    `z value`= round(`z value`, 3),
    `p value`= signif(`Pr(>|z|)`, 3)
  )
print(Table_S7, row.names = FALSE)

# VISUALISE PREDICTIONS
use_mod <- MuMIn::get.models(ms, 1)[[1]]
model7_elev_night_rain <- use_mod

re1 <- try(ranef(model7_elev_night_rain)[[1]], silent = TRUE)
if (!inherits(re1, "try-error")) {
  if (is.data.frame(re1)) { re_vec <- re1[,1]; names(re_vec) <- rownames(re1) } else { re_vec <- re1 }
  site0 <- names(re_vec)[which.min(abs(as.numeric(re_vec)))]
} else {
  site0 <- levels(df_abund$site_ID)[1]
}
site0 <- factor(site0, levels = levels(df_abund$site_ID))

cat_levels <- levels(df_abund$bib_abundance_cat)

pred_one_cat <- function(newdat, k) {
  newdat$bib_abundance_cat <- factor(cat_levels[k], levels = cat_levels, ordered = TRUE)
  as.numeric(predict(model7_elev_night_rain, newdata = newdat, type = "prob"))
}

## ELEVATION curve (others at mean: z = 0; Night = 1)
pred_dat <- data.frame(
  elevation = seq(min(df_long$elevation, na.rm = TRUE),
                  max(df_long$elevation, na.rm = TRUE), length.out = 100),
  next_day_precip_std = 0,
  night = 1
)
pred_dat$elevation_std <- (pred_dat$elevation - mean(df_long$elevation, na.rm=TRUE)) /
  sd(df_long$elevation, na.rm=TRUE)
pred_dat$site_ID <- site0

p1 <- pred_one_cat(pred_dat, 1)
p2 <- pred_one_cat(pred_dat, 2)
p3 <- pred_one_cat(pred_dat, 3)
p4 <- pred_one_cat(pred_dat, 4)

elevation_preds <- data.frame(
  Elevation   = pred_dat$elevation,
  Abundance_1 = p1,
  Abundance_2 = p2,
  Abundance_3 = p3,
  Abundance_4 = p4
)

## NIGHT vs DAY bars (hold others at z = 0)
pred_dat <- data.frame(
  elevation_std       = 0,
  next_day_precip_std = 0,
  night               = c(0, 1),
  site_ID             = site0
)

p1 <- pred_one_cat(pred_dat, 1)
p2 <- pred_one_cat(pred_dat, 2)
p3 <- pred_one_cat(pred_dat, 3)
p4 <- pred_one_cat(pred_dat, 4)

night_preds <- data.frame(
  Night        = pred_dat$night,
  Abundance_1  = p1,
  Abundance_2  = p2,
  Abundance_3  = p3,
  Abundance_4  = p4
)

## RAIN curve (others at z = 0; Night = 1) 
pred_dat <- data.frame(
  elevation_std       = 0,
  next_day_precip     = seq(min(df_long$next_day_precip, na.rm=TRUE),
                            max(df_long$next_day_precip, na.rm=TRUE), length.out = 100),
  night               = 1,
  site_ID             = site0
)
mean_rain <- mean(df_long$next_day_precip, na.rm=TRUE)
sd_rain   <- sd(df_long$next_day_precip, na.rm=TRUE)
pred_dat$next_day_precip_std <- (pred_dat$next_day_precip - mean_rain) / sd_rain

p1 <- pred_one_cat(pred_dat, 1)
p2 <- pred_one_cat(pred_dat, 2)
p3 <- pred_one_cat(pred_dat, 3)
p4 <- pred_one_cat(pred_dat, 4)

rain_preds <- data.frame(
  Rain         = pred_dat$next_day_precip,
  Abundance_1  = p1,
  Abundance_2  = p2,
  Abundance_3  = p3,
  Abundance_4  = p4
)

# Check column names just to confirm what's in elevation_preds
colnames(elevation_preds)

elevation_long <- elevation_preds %>%
  pivot_longer(
    cols = starts_with("Abundance_"),
    names_to = "Abundance",
    values_to = "Probability"
  )

unique(elevation_long$Abundance)

elevation_long <- elevation_long %>%
  mutate(Abundance = case_when(
    Abundance == "Abundance_1" ~ "1 (1–5 individuals)",
    Abundance == "Abundance_2" ~ "2 (6–20 individuals)",
    Abundance == "Abundance_3" ~ "3 (21–50 individuals)",
    Abundance == "Abundance_4" ~ "4 (51–100 individuals)"
  ))

p_elev <- ggplot(elevation_long, aes(x = Elevation, y = Probability, color = Abundance)) +
  geom_line(size = 1.2) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Elevation (m)", y = "Predicted Probability", tag = "B") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")

rain_long <- rain_preds %>%
  pivot_longer(
    cols = starts_with("Abundance_"),
    names_to = "Abundance",
    values_to = "Probability"
  )

rain_long <- rain_long %>%
  mutate(Abundance = case_when(
    Abundance == "Abundance_1" ~ "1 (1–5 individuals)",
    Abundance == "Abundance_2" ~ "2 (6–20 individuals)",
    Abundance == "Abundance_3" ~ "3 (21–50 individuals)",
    Abundance == "Abundance_4" ~ "4 (51–100 individuals)"
  ))

p_rain <- ggplot(rain_long, aes(x = Rain, y = Probability, color = Abundance)) +
  geom_line(size = 1.2) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Next-day Rainfall (mm)", y = "Predicted Probability", tag = "A") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")  # Hide to avoid duplication

p_rain

library(patchwork)

combined_plot <- p_elev + plot_layout(ncol = 1)
combined_plot

ggsave("Figure6_abundance_predictors_combined.png", combined_plot, width = 15, height = 5, dpi = 300)


night_long <- night_preds %>%
  pivot_longer(
    cols = starts_with("Abundance_"),
    names_to = "Abundance",
    values_to = "Probability"
  )

night_long <- night_long %>%
  mutate(Abundance = case_when(
    Abundance == "Abundance_1" ~ "1 (1–5 individuals)",
    Abundance == "Abundance_2" ~ "2 (6–20 individuals)",
    Abundance == "Abundance_3" ~ "3 (21–50 individuals)",
    Abundance == "Abundance_4" ~ "4 (51–100 individuals)"
  ))

night_long$Night <- factor(night_long$Night, levels = c(0, 1), labels = c("Day", "Night"))

#bar graph 
p_night <- ggplot(night_long, aes(x = Night, y = Probability, fill = Abundance)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "", y = "Predicted Probability") +
  theme_classic(base_size = 14)

p_night

ggsave("Figure6_abundance_predictors_night.png", p_night, width = 6, height = 5, dpi = 300)

combined_plot_1col <- (p_elev / p_night) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

combined_plot_1col

ggsave("Figure_3.png", combined_plot_1col, width = 8, height = 10, dpi = 400)



