#GRAPHING OREI PLFAs
  #This script creates PCA plots on PLFA data from OREI 2021 samples and tests for differences between   groups.
  #Hannah Rodgers 3/23/23

library(readxl)
library(vegan)
library(missMDA)
library(patchwork)
library(ggrepel)
library(tidyverse)

#load data
OREI_PLFA <- read_excel("data/OREI_PLFA_2021.xlsx")

#make pca object on raw biomarkers
PCA <- prcomp(OREI_PLFA[,c(24:46)], center = TRUE, scale. = TRUE)

#put PC1 and PC2 into dataframe
OREI_PLFA[, c('PC1', 'PC2')] = PCA$x[, 1:2]

#use envfit to correlate other PLFA groupings to ordination
vars <- select(OREI_PLFA, Fungi, Bacteria, AMF)

en <- vegan::envfit(PCA, vars, 
                    permutations = 10000, na.rm = TRUE, strata = NULL)

# use results from envfit to create arrows (continuous variables)
arrows_cont = as.data.frame(scores(en, "vectors"))

# use results from envfit to create arrows (categorical variables)
#en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#rescale arrows to fit nicely within scatterplot of our data
arrows_cont$var = rownames(arrows_cont)

mult = max(abs(OREI_PLFA[, c('PC1', 'PC2')])) / max(abs(arrows_cont[, 1:2])) / 2
arrows_cont[, 1:2] = arrows_cont[, 1:2] * mult

#use PERMANOVA (adonis testing) to test for sig differences between groups
adonis2(OREI_PLFA[,c(22:44)] ~ Compost_Rate * Compost_Rate_Year * Cover_Crop * Compost_Application_Year * Crop_Phase, data = OREI_PLFA)

#test for sig differences in dispersion (distance from centroid) between groups
disp <- betadisper(dist(OREI_PLFA[,c(22:44)]),
                   #here, chose the grouping variable or use paste to group by 2 variables
                   OREI_PLFA$Cover_Crop, 
                   type = 'centroid')

anova(disp)
boxplot(disp)

#reorder compost rate
OREI_PLFA$Compost_Rate <- factor(OREI_PLFA$Compost_Rate, 
        levels = c("Control", "Inorganic Fertilizer", "Low", "Medium", "High"))

#### PLOT IT! Use OREI_PLFAs for scatterplot, and arrows_cont for arrows ####

#COMPOST RATE AND YEARS SINCE APPLICATION
(COMPOST <- OREI_PLFA %>% 
  ggplot(aes(x = PC1, y = PC2, 
             linetype = Years_Since_Application, 
             shape = Years_Since_Application, 
             color = Compost_Rate)) +
  
  geom_point(size = 1.8) +
  scale_color_manual (values = c("#8C510A","#DFC27D","#5BB300","#00B4EF","#CF78FF")) +
  
  stat_ellipse(level = 0.95, lwd = 0.8) +
  
  geom_segment(data = arrows_cont, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2), inherit.aes = FALSE,
               arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = arrows_cont, 
                   aes(PC1 * 1.1, PC2 * 1.1, label = var), 
                   inherit.aes = FALSE, 
                   nudge_x = -0.5, 
                   box.padding = 0.1) +
  
  theme_bw(base_size = 13) +
  labs(colour = "Compost Rate", 
       linetype = "Years Since Application", 
       shape = "Years Since Application") +
   
   guides (color = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE),
           linetype = guide_legend(override.aes = list(lwd = 0.5))) +
   
   theme(legend.position = "bottom", legend.box = "vertical", 
         legend.background = element_rect(colour = "black", fill = "white", linetype = "solid")))


#COVER CROPS AND CROP PHASE
(COVER <- OREI_PLFA %>% 
  ggplot(aes(x = PC1, y = PC2, 
             linetype = Cover_Crop, 
             shape = Cover_Crop, 
             color = Crop_Phase)) +
  geom_point(size = 1.8) +
  stat_ellipse(level = 0.95, lwd = 1) +
  geom_segment(data = arrows_cont, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2), inherit.aes = FALSE,
               arrow = arrow(length = unit(0.02, "npc"))) +
  
  geom_label_repel(data = arrows_cont, 
                   aes(PC1 * 1.1, PC2 * 1.1, label = var), 
                   inherit.aes = FALSE, 
                   nudge_x = -0.5, 
                   box.padding = 0.1) +
  
  theme_bw(base_size = 13) +
  labs(colour = "Crop Phase", 
       linetype = "Cover Crops", 
       shape = "Cover Crops") +
  guides (linetype = guide_legend(override.aes = list(lwd = 0.5))) +
    
    theme(legend.position = "bottom", legend.box = "vertical", 
          legend.background = element_rect(colour = "black", fill = "white", linetype = "solid")))

#COMBINE GRAPHS WITH PATCHWORK
#save as 1023*665
COMPOST + COVER + 
  plot_annotation(tag_levels = 'a')

####OTHER THINGS####
#find which variables contribute most to the PCA
  res.var <- get_pca_var(PLFA_pca)
  res.var$contrib
  
  
#STRESS RATIOS FROM VAN DIEPEN ET AL, 2010
    #only using ratio of 17:0cyclo to 16:1.w7c, because I don't have any 19:0cyclo

OREI_PLFA$stress <- OREI_PLFA$X17.0.cyclo / OREI_PLFA$X16.1.w7c

OREI_PLFA$MRR <- OREI_PLFA$X17.0.cyclo / OREI_PLFA$X16.0

OREI_PLFA$BRR <- OREI_PLFA$X16.0.10.methyl / OREI_PLFA$X16.0

#linear regression of ratios and compost
OREI_PLFA$Compost_2015 <- as.numeric(OREI_PLFA$Compost_2015)
OREI_PLFA$Compost_2020 <- as.numeric(OREI_PLFA$Compost_2020)
OREI_PLFA <- OREI_PLFA %>% 
  filter(Compost_Rate != "Inorganic Fertilizer")

mod <- lm(BRR ~ as.numeric(Compost_2020), data = OREI_PLFA)
summary(mod)
summary(mod)$coefficients[2,1]
plot(stress ~ Compost_2020, data = OREI_PLFA)


#create new column with total compost

                            
OREI_PLFA$allcompost <- rowSums(OREI_PLFA[,c("Compost_2015", "Compost_2020")], na.rm=TRUE)

#duplicate controls! change one set to 2020, other set to 2016
OREI_Xs <- OREI_PLFA[which(OREI_PLFA$Years_Since_Application == 'none'), ] %>% 
  mutate(Years_Since_Application = '1 year')

OREI_duped <- OREI_PLFA %>% 
  mutate(Years_Since_Application = ifelse(Years_Since_Application == 'none', '5-6 years',
                               ifelse(Years_Since_Application == '5-6 years', '5-6 years',
                                      ifelse(Years_Since_Application == '1 year', '1 year', NA)))) %>% 
  bind_rows(OREI_Xs) %>% 
  filter(Compost_Rate != "Inorganic Fertilizer")

library(ggpubr)

#PLOT IT!
ggplot(OREI_duped, aes(x= allcompost, y = BRR, color = Years_Since_Application)) +
  geom_point() +
  theme_bw(base_size = 16) +
  geom_smooth(method = "lm") +
  labs(x = "Compost Rate", y = "Gram - Stress and Carbon Limitation \n (ratio of 17:0cyclo to 16.1.w7c)",
       color= "Years Since \nCompost Application") +
stat_cor() 


                        