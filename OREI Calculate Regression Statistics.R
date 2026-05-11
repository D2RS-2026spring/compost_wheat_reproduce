#This R script is used to clean up and transform soil data to perform regression analysis.
#author: Hannah Rodgers

## Load Packages
library(MASS)
#library(car)
library(readxl)
#library(ggpubr)
library(tidyverse)

#### IMPORT AND TIDY DATA ####

#run this script for each group: 2021 wheat, 2021 fallow, 2022 fallow, 2022 wheat

#select the variable we're interested in
vars <- c('Sampling_Year', 'Total_Compost', 'Years_Since_Application', 'Crop_Phase',
          'POXC', 'PMN', 'DOC', 'DON', 'PMC', 'Protein', 'TN', "SOC")

## Import Data, remove IF and compost not_yet_applied, keep only variables in vars
OREI_2021 <- read_excel("data/OREI_soil_2021.xlsx") %>% 
  filter(Compost_Rate != "not_yet_applied", 
         Compost_Rate != "Inorganic Fertilizer") %>% 
 select(all_of(vars))

#linear regression (to check that its as expected)
  mod <- lm(SOC ~ Total_Compost, data = OREI_2021)
  summary(mod)
  plot(SOC ~ Total_Compost, data = OREI_2021)

#separate by phase and years_since_application. Drop columns w only NAs (for fallow)
  #Ignore fallow 2021 since not all compost had been applied. Run one year and phase at a time.
w_22_1y <- subset(OREI_2022, Crop_Phase == "wheat" & 
                  Years_Since_Application != "5-6 years")

w_22_5y <- subset(OREI_2022, Crop_Phase == "wheat" & 
                  Years_Since_Application != "1 year")

f_22_5y <- subset(OREI_2022, Crop_Phase == "fallow" & 
              #    Years_Since_Application != "1-2 years")%>%
        #  select(where(function(x) any(!is.na(x))))

#### CALCULATE REGRESSION STATISTICS ####
  #run this script on each group separately (save it as data)
  #we'll check regression assumptions and adjust p-values later
data <- f_21_5y

#this function saves R2, slope, st error, p-value, and SLOPE/MEAN * 50 (% change per 50 Mg compost)
lin_reg <- function(x) {
  mod <- lm(x ~ Total_Compost, data = data)
  
  #save slope, r2, and p-value of regression
  my_output <- c(summary(mod)$coefficients[2,1],
                 summary(mod)$r.squared,
                 summary(mod)$coefficients[2,4],
                 
  #save std error (/ mean, * 50 * 100 to match % change data)
                 summary(mod)$coefficients[2,2]/mean(x, na.rm = TRUE) * 50 * 100)
  
  #save % change: divide slope by mean, * 50 * 100 to get % change for 50 Mg/ha compost
  my_output[5] <- (my_output[1]/mean(x, na.rm = TRUE)) * 50 * 100
  return(my_output)}

#RUN THE FUNCTION on all variables
reg_stats <- as.data.frame(lapply(data[,c(5:ncol(data))], lin_reg))

#clean up reg_stats
reg_stats <- reg_stats %>% 
  t() %>% 
  as.data.frame()

colnames(reg_stats) <- c("Slope", "R_Squared", "P_Value", "St_Error", "Percent_Change")
reg_stats$variable <- as.factor(rownames(reg_stats))

#which data is this?
head(data)

#add that metadata to regstats_phase_year
reg_stats_f_21_5y <- reg_stats %>% 
  mutate(Years_Since_Compost_Application = "5-6 years",
         Crop_Phase = "Fallow",
         Sampling_Year = "2021")

#STOP HERE and run other group (1y or 5y)

#Bind 1yr and 5yr together
reg_stats <- bind_rows(reg_stats_w_22_1y, reg_stats_w_22_5y)

#add stars from p values to that df
reg_stats$stars <- ifelse(reg_stats$P_Value <= 0.001, "***",
                              ifelse(reg_stats$P_Value <= 0.01, "**",
                                     ifelse(reg_stats$P_Value <= 0.05, "*", "")))

#Save
write.csv(reg_stats, "data/reg_stats/reg_stats_fallow_2021.csv")

#### CHECK REGRESSION ASSUMPTIONS ####
#Run each group through this script separately. If data needs to be transformed, manually update p-values in reg_states dataframe
data <- f_21_5y

#CHECK NORMALITY OF RESIDUALS
#Create an empty list
Shapiro.pvals <- list()

#This loop runs a linear model on compost ~ each variable, 
#evaluates normality with the Shapiro test, and saves the p value.
for (i in names(data[,c(5:ncol(data))])) {
  mod <- lm(get(i) ~ Total_Compost, data = data)
  Shapiro.pvals[[i]] <- (shapiro.test(mod$residuals))$p.value }

#print variables with p < 0.05
as.data.frame(t(as.data.frame(Shapiro.pvals))) %>% 
  filter(V1 < 0.05)

## for non-normal variables, box-cox transform, re-check normality, then calculate updated p-value #
v <- (data$DON)^-0.2
mod <- lm(v ~ Total_Compost, data = data)

boxcox(mod)

hist(resid(mod))
shapiro.test(resid(mod))

summary(mod)




## CHECK LINEARITY WITH F-TEST

#This function creates a linear and quadratic model for each variable, 
#tests whether the two models differ using ANOVA, and saves the p.value
F_test <- function(x) {
  mod <- lm(x ~ Total_Compost, data = data)
  reduced<-lm(x ~ Total_Compost, data = data)
  full<-lm(x ~ Total_Compost, data = data)
  return(anova(reduced, full)$"Pr(>F)") }

#run this function on each variable and save the p.values in a list
F.test <- list()
for (i in colnames(data[,c(12:ncol(data))])) {
  F.test[[i]] <- F_test(data[[i]]) }

#print variables with p < 0.05
as.data.frame(t(as.data.frame(F.test))) %>% 
  filter(V1 < 0.05)

## CHECK CONSTANCY OF RESIDUALS WITH THE LEVENE TEST ##
#run levene test and save p-values
levene <- list()
for (i in colnames(data[,c(5,12:ncol(data))])) {
  result <- leveneTest((data[[i]]) ~ as.factor(data$Total_Compost))
  levene[[i]] <- result$`Pr(>F)`[1] }

#print variables with p < 0.05
as.data.frame(t(as.data.frame(levene))) %>% 
  filter(V1 < 0.05)
