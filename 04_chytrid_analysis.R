
##########---------------CHYTRID-----------------########

#### CHYTRID PREVELANCE AND INTENSITY 
library(readxl)
library(dplyr)
library(ggplot2)
library(ggthemes)

chytrid_data <- read_excel("Bd_summary.xlsx")

# Clean column names 
names(chytrid_data) <- trimws(names(chytrid_data))

chytrid_data <- chytrid_data %>%
  mutate(
    Zoospore_Eq = as.numeric(Zoospore_Eq),
    Bd = as.numeric(Bd)
  )

# Summarise per site
chytrid_summary <- chytrid_data %>%
  group_by(Site) %>%
  summarise(
    Mean_Zoospore_Eq = if (all(is.na(Zoospore_Eq))) NA_real_ else mean(Zoospore_Eq, na.rm = TRUE),
    Max_Zoospore_Eq = if (all(is.na(Zoospore_Eq))) NA_real_ else max(Zoospore_Eq, na.rm = TRUE),
    Min_Zoospore_Eq = if (all(is.na(Zoospore_Eq))) NA_real_ else min(Zoospore_Eq, na.rm = TRUE),
    Bd_Prevalence = if (all(is.na(Bd))) NA_real_ else mean(Bd, na.rm = TRUE) * 100,
    Elevation = mean(ELEV, na.rm = TRUE),
    Count = n()
  )

chytrid_summary <- chytrid_summary %>%
  mutate(Log_Zoospore_Eq = log1p(Mean_Zoospore_Eq))

chytrid_summary <- chytrid_summary %>%
  mutate(Elevation_Category = cut(Elevation, breaks=c(0, 500, 999, 1800), labels=c("Low (0-500 m asl)", "Mid (501-999 m asl)", "High (1000-1700 m asl)")))


write.csv(chytrid_summary, "Bd_infection_summary_new_elevation_groups.csv", row.names = FALSE)
