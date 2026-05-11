## Load Packages
library(readxl)
library(lme4)
library(emmeans)
library(multcomp)
library(tidyverse)
library(car)

#### IMPORT AND TIDY DATA ####

## Import Data (first OREI_2021, then OREI_2022)
OREI_2021 <- read_excel("data/OREI_soil_2021.xlsx")
OREI_2022 <- read_excel("data/OREI_soil_2022.xlsx")

OREI_Enzymes <- read_excel("data/OREI_enzymes.xlsx")

#bind together, remove plots inorganic_fertilizer or compost not_yet_applied
OREI <- bind_rows(OREI_2021, OREI_2022) %>% 
  mutate(Total_Compost = as.factor(Total_Compost)) %>% 
  filter(Compost_Rate != "not_yet_applied", 
         Compost_Rate != "Inorganic Fertilizer")

#### CALCUALTE MEANS AND SD ####
OREI %>% 
  filter(Compost_2015 == 45 | Compost_2020 == 50 | Compost_2020 == 0) %>% 
  group_by (Compost_Rate_Year) %>% 
  ggpubr::get_summary_stats(DON, type = "mean_sd")
  
#### P VALUES FOR OVERALL EFFECTS AND INTERACTIONS USING LMER####
  #for first table, use wheat only
   OREI_wheat <- subset(OREI, Crop_Phase == "wheat") 

#test for differences and interactions
OREI$var <- (OREI$InorganicC)

mod <- lmer(var ~ Cover_Crops * Compost_Rate_Year +
               (1|Sampling_Year), 
             data=OREI, REML=T)

#test for normality
hist(resid(mod))
shapiro.test(resid(mod))

#if not normal, log transform (must use aov to get boxcox)
      mod2 <- lm((var+1) ~ Cover_Crops * Compost_Rate_Year * Sampling_Year, 
             data=OREI)
      MASS::boxcox(mod2)

#get p-values for comparing means
cld(emmeans(mod, ~ Cover_Crops), 
              Letters = LETTERS, adjust = "none")

#get p-values for main effects and interactions
Anova(mod)
