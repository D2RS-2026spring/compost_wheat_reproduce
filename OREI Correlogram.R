#GRAPHING CORRELOGRAM
  #This script creates a correlation matrix to look at correlations between soil microbiology and other soil properties.
  #Hannah Rodgers 3/23/23

library(readxl)
library(tidyverse)
library(ggcorrplot)

#load all data
OREI_2021 <- read_xlsx("data/OREI_soil_2021.xlsx")
OREI_2022 <- read_xlsx("data/OREI_soil_2022.xlsx")
OREI_PLFA <- read_xlsx("data/OREI_PLFA_2021.xlsx")
OREI_enzymes <- read_xlsx("data/OREI_enzymes.xlsx")

#merge dataframes
OREI <- bind_rows(OREI_2021, OREI_2022)
OREI <- left_join(OREI_enzymes, OREI, by = 'Sample_ID_Year')
OREI <- left_join(OREI, OREI_PLFA, by = 'Sample_ID_Year')

#use geometric mean to calculate C_enzymes, N_enzymes
library(compositions)
C_enzymes <- select(OREI, AG, BG, BX, CBH)
N_enzymes <- select(OREI, LAP, NAG)

OREI$Carbon_Enzymes <- geometricmeanRow(C_enzymes)
OREI$Nitrogen_Enzymes <- geometricmeanRow(N_enzymes)

#select data to correlate
vars <- dplyr::select(OREI,
                      Carbon_Enzymes, Nitrogen_Enzymes, #enzymes
                      SOC, TN, POXC, PMC, DOC, DON, PMN, Protein, #soil health
                      H2O, #soil chem
                      Fungi, Bacteria) #PLFAs

#rename
vars <- rename(vars,
               'C Enzymes' = Carbon_Enzymes,
               'N Enzymes' = Nitrogen_Enzymes,
               "Total N" = TN,
              'SOC' = SOC,
               'POXC' = POXC,
               'Soil Moisture' = H2O,
               'DOC' = DOC, 
                'DON'  = DON,
               'PMN' = PMN,
               'Soil Protein' = Protein,
               'Fungal Biomass' = Fungi,
               'Bacterial Biomass' = Bacteria)

#make correlation matrix and prune it
corr <- cor(vars, use = 'pairwise.complete.obs')

corr_pruned <- corr %>% 
  as.data.frame() %>% 
  select("C Enzymes", "N Enzymes", 'Fungal Biomass', 'Bacterial Biomass')

corr_pruned <- corr_pruned[c(-1,-2),]

corr_pruned$`Fungal Biomass`[10:11]=0
corr_pruned$`Bacterial Biomass`[10:11]=0

#make p value matrix
p.mat <- cor_pmat(vars)

p.mat_pruned <- p.mat %>% 
  as.data.frame() %>% 
  select("C Enzymes", "N Enzymes", 'Fungal Biomass', 'Bacterial Biomass') %>% 
  as.matrix()

p.mat_pruned <- p.mat_pruned[c(-1,-2),]

#create correlogram
ggcorrplot(corr_pruned, p.mat = p.mat_pruned, insig = "pch", outline.color = "darkslategrey",
           legend.title = "Strength of \nCorrelation (R)",
           ggtheme = ggplot2::theme_minimal(base_size = 16, element_text(color = "black")))

####OTHER

#calculate new parameters: C/PLFA, C/fungi, C/bacteria
OREI$C_Enzymes_to_all_PLFA <- OREI$Carbon_Enzymes/OREI$Total_micro
OREI$C_Enzymes_to_fungi <- OREI$Carbon_Enzymes/OREI$Fungi
OREI$C_Enzymes_to_bacteria <- OREI$Carbon_Enzymes/OREI$Bacteria

#analyze new parameters
mod <- lm(C_Enzymes_to_bacteria ~ Compost_2020.x, data = OREI)
summary(mod)
plot(C_Enzymes_to_all_PLFA ~ Compost_2015.x, data = OREI)

#plot new parameters
OREI %>% 
  ggplot(aes(x = C_Enzymes_to_all_PLFA, 
             y = Total_Compost.x,
  color = Years_Since_Application))+
  geom_smooth(method = 'lm', size = 3, se = F) +
  stat_regline_equation(aes(label =  paste(..rr.label..)), size = 7)
