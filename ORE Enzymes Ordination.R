#GRAPHING OREI ENZYME ACTIVITIES
  #This script creates PCA plots on Enzyme data from OREI 2021 & 2022 samples and tests for differences   between   groups.
  #Hannah Rodgers 3/23/23

library(readxl)
library(vegan)
library(missMDA)
library(patchwork)
library(ggrepel)
library(tidyverse)

#load data
OREI_enzymes <- read_excel("Data/OREI_enzymes.xlsx")

#impute NAs, then merge imputed data back into dataframe
imputed <- imputePCA(OREI_enzymes[,c(12:19)])
imputed <- as.data.frame(imputed$completeObs)

OREI_enzymes <- OREI_enzymes %>% 
  select(-c(AG:SUL)) %>% 
  bind_cols(imputed)

#create PCA object
PCA     <- prcomp(OREI_enzymes[,c(12:19)], center = TRUE, scale = TRUE)

#put PC1 and PC2 into dataframe
OREI_enzymes  [, c('PC1', 'PC2')] = PCA$x[, 1:2]

# save variable loadings in a separate dataframe called rot
rot = as.data.frame(PCA$rotation[, 1:2])
rot$var = rownames(PCA$rotation)

# rescale the loadings to fit nicely within the scatterplot of our data
mult = max(abs(OREI_enzymes[, c('PC1', 'PC2')])) / max(abs(rot[, 1:2])) / 2
rot[, 1:2] = rot[, 1:2] * mult

#use PERMANOVA (adonis testing) to test for significant differences between groups
adonis2(OREI_enzymes[,c(12:19)] ~ Compost_Rate_Year * Cover_Crops, data = OREI_enzymes)

#test for sig differences in dispersion (distance from centroid) between groups
OREI_enzymes[is.na(OREI_enzymes)] <- 0
OREI_enzymes_w <- filter(OREI_enzymes, Crop_Phase == "Fallow")

disp <- betadisper(dist(OREI_enzymes[,c(12:19)]), 
                   #here, chose the grouping variable or use paste to group by 2 variables
                  OREI_enzymes$Cover_Crops,
                  type = "centroid")

anova(disp)
boxplot(disp)
plot(disp)

#reorder the levels of compost rates and save sampling year as a factor
OREI_enzymes$Compost_Rate<- factor(OREI_enzymes$Compost_Rate, 
        levels = c("None", "Inorganic Fertilizer", "Low", "Medium", "High" ))

OREI_enzymes$Sampling_Year <- as.character (OREI_enzymes$Sampling_Year)

#### PLOT IT! Use OREI_enzymes for scatterplot, and rot for arrows ####

#COMPOST RATE AND YEARS SINCE APPLICATION
(COMPOST <- OREI_enzymes %>% 
ggplot(aes(x = PC1, y = PC2, 
           color = Compost_Rate,
           linetype = Years_Since_Application, 
           shape = Years_Since_Application)) +
  
  geom_point(size = 1.8) +
  scale_color_manual (values = c("#8C510A","#DFC27D","#5BB300","#00B4EF","#CF78FF")) +
  
  stat_ellipse(level = 0.95, lwd = 1) +
  
  geom_segment(data = rot, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2), inherit.aes = FALSE, 
               arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = rot, 
             aes(PC1 * 1.1, PC2 * 1.1, label = var), 
             inherit.aes = FALSE, 
             nudge_x = 0.5, 
             box.padding = 0.1) +
  
  theme_bw(base_size = 13)+
  labs(colour = "Compost Rate", 
       linetype = "Years Since Application", 
       shape = "Years Since Application") + 

  guides (color = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE),
          linetype = guide_legend(override.aes = list(lwd = 0.5))) +
  
  theme(legend.position = "bottom", legend.box = "vertical", 
        legend.background = element_rect(colour = "black", fill = "white", linetype = "solid")))

#CROPPING PHASE AND SAMPLING YEAR
(PHASE <- OREI_enzymes %>% 
  ggplot(aes(x = PC1, y = PC2, 
             linetype = Sampling_Year, 
             shape = Sampling_Year, 
             color = Crop_Phase)) +
  
  geom_point() +
  scale_shape_manual(values=c(17,16))+
  
  stat_ellipse(level = 0.95, lwd = 0.9) +
  geom_segment(data = rot, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2), inherit.aes = FALSE, 
               arrow = arrow(length = unit(0.02, "npc"))) +
  
  geom_label_repel(data = rot, 
                   aes(PC1 * 1.1, PC2 * 1.1, label = var), 
                   inherit.aes = FALSE, 
                   nudge_x = 0.5, 
                   box.padding = 0.1) +
  
  theme_bw(base_size = 13)+
  labs(colour = "Cropping Phase", 
       linetype = "Sampling Year", 
       shape = "Sampling Year") +
  
  theme(legend.position = "bottom", legend.box = "vertical", 
        legend.background = element_rect(colour = "black", fill = "white", linetype = "solid")))

#COMBINE GRAPHS WITH PATCHWORK
#save as 1023*665
COMPOST + PHASE +
  plot_annotation(tag_levels = 'a')

#find which variables contribute most to the PCs
library(factoextra)
res.var <- get_pca_var(PCA)
res.var$contrib

#### C v N v P RATIO GRAPHS ####
OREI <- read_xlsx("CHECKED data/OREI_Enzymes.xlsx")

#order levels of compost rate and save sampling year as a factor
OREI$Compost_Rate<- factor(OREI$Compost_Rate, 
              levels = c("None", "Inorganic Fertilizer", "Low", "Medium", "High" ))

OREI$Sampling_Year<- as.factor(OREI$Sampling_Year)

# Create C:N and C:P proportions
OREI <- OREI %>% 
  mutate(CN = BG / (BG + NAG + LAP), 
         CP = BG / (PHOS + BG))

#test significance of new variables (CN and CP) with anova
OREI_2021 <- filter(OREI, Sampling_Year == "2021")
OREI_2022 <- filter(OREI, Sampling_Year == "2022")

summary(aov(OREI_2021$CN ~ Compost_Rate_Year * Cover_Crops * Rotation, data = OREI_2021))

#GRAPH: COMPOST RATE AND APPLICATION YEAR
OREI %>% 
  ggplot(aes(CP, CN, 
             color = Compost_Rate, 
             linetype = Years_Since_Application,
             shape = Years_Since_Application)) +
  stat_ellipse(level = 0.95) + 
  geom_point() +
  scale_color_manual (values = c("#8C510A","#DFC27D","#5BB300","#00B4EF","#CF78FF")) +
  xlim(0, 1) + ylim(0, 1) +
  # hline and vline add the quadrants
  geom_hline(yintercept = 0.5, alpha = 0.5) + 
  geom_vline(xintercept = 0.5, alpha = 0.5) +
  annotate(geom = "text", x = 0.1, y = 1, label = "P Limited") +
  annotate(geom = "text", x = 0.1, y = 0.1, label = "N/P Limited") + 
  annotate(geom = "text", x = .9, y = 1, label = "C Limited") +
  annotate(geom = "text", x = .9, y = 0.1, label = "N Limited") +
  ylab(label = "C:N") + xlab(label = "C:P") +
  ggtitle(label = "C vs. N vs. P Acquisition") +
  # labs(subtitle = "Post-Treatment") + 
  # All the code below here is just personal preference for how plot looks
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank()) +
  theme(plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
        legend.text = element_text(size =11, colour = "grey30"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(colour = "Compost Rate", 
       linetype = "Years Since Application", 
       shape = "Years Since Application") + 
  guides (linetype = guide_legend(override.aes = list(lwd = 0.5)))


#SAMPLING YEAR AND ROTATION
OREI %>% 
  ggplot(aes(CP, CN, 
             color = Rotation, 
             linetype = Sampling_Year,
             shape = Sampling_Year)) +
  stat_ellipse(level = 0.95) + 
  geom_point() +
  xlim(0, 1) + ylim(0, 1) +
  # hline and vline add the quadrants
  geom_hline(yintercept = 0.5, alpha = 0.5) + 
  geom_vline(xintercept = 0.5, alpha = 0.5) +
  annotate(geom = "text", x = 0.1, y = 1, label = "P Limited") +
  annotate(geom = "text", x = 0.1, y = 0.1, label = "N/P Limited") + 
  annotate(geom = "text", x = .9, y = 1, label = "C Limited") +
  annotate(geom = "text", x = .9, y = 0.1, label = "N Limited") +
  ylab(label = "C:N") + xlab(label = "C:P") +
  ggtitle(label = "C vs. N vs. P Acquisition") +
  # labs(subtitle = "Post-Treatment") + 
  # All the code below here is just personal preference for how plot looks
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_blank()) +
  theme(plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"),
        legend.text = element_text(size =11, colour = "grey30"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(colour = "Crop Phase", 
       linetype = "Sampling Year", 
       shape = "Sampling Year") + 
  guides (linetype = guide_legend(override.aes = list(lwd = 0.5)))
