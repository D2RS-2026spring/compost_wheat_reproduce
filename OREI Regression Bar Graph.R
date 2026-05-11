library(readxl)
library(ggpubr)
library(patchwork)
library(tidyverse)

#### DATA INPUT ####
#read in files
reg_stats_f_21 <- read_csv("data/reg_stats/reg_stats_fallow_2021.csv")
reg_stats_f_22 <- read_csv("data/reg_stats/reg_stats_fallow_2022.csv")
reg_stats_w_21 <- read_csv( "data/reg_stats/reg_stats_wheat_2021.csv")
reg_stats_w_22 <- read_csv( "data/reg_stats/reg_stats_wheat_2022.csv")

OREI_2021 <- read_excel("data/OREI_soil_2021.xlsx")

#### SCATTERPLOTS ####

OREI_2021 <- OREI_2021 %>% 
  filter(Compost_Rate != "not_yet_applied", 
         Compost_Rate != "Inorganic Fertilizer", 
         Crop_Phase == "wheat")

#duplicate controls!
OREI_Xs <- OREI_2021[which(OREI_2021$Years_Since_Application == 'none'), ] %>% 
  mutate(Years_Since_Application = '1 year')

OREI_duped <- OREI_2021 %>% 
  mutate(Years_Since_Application = 
           ifelse(Years_Since_Application == 'none', '5-6 years',
           ifelse(Years_Since_Application == '5-6 years', '5-6 years',
           ifelse(Years_Since_Application == '1 year', '1 year', NA)))) %>% 
  bind_rows(OREI_Xs)

## make regression scatterplot ##
(scatter <- OREI_duped %>% 
  ggplot(aes(x = Total_Compost, 
             y = DON, 
             color = Years_Since_Application))+
  geom_point(size = 3) +
  geom_smooth(method = 'lm', size = 3, se = F) +
  stat_regline_equation(aes(label =  paste(..rr.label..)), size = 6) +
  labs(x = "Compost (Mg/ha)", 
       y = "DOC (mg/kg)",
       color= "Years Since \nCompost Application") +
  theme_bw(base_size = 15) + theme(legend.position = "none"))

#### PERCENT CHANGE GRAPH ####

#reorder bars
#reg_stats$variable_name <- factor(reg_stats_om$variable_name, 
                 # levels = c("Dissolved C", "Dissolved N", 
                    #         "Mineralizable C", "Mineralizable N", "Oxidizable C", 
                    #         "Protein", "Total Organic C", "Total N", "Aggregate Stability"))

#PLOT IT
(p3 <- reg_stats_f_22 %>% 
  
ggplot( aes(y = Percent_Change, 
            x = variable, 
            fill = Years_Since_Compost_Application)) +
  geom_bar(stat = "identity", 
           position = "dodge", size = 3) +
  
  geom_errorbar(aes(ymin = Percent_Change - St_Error, 
                    ymax = Percent_Change + St_Error), 
                position = position_dodge(0.9), width = 0.2) + 
  
  labs(title = '2022 Fallow',
       fill = "Years Since\nCompost Application") +
  
  theme_bw(base_size = 15) +
  theme(axis.title = element_blank()) +
  #theme(legend.position = "none") +
  
  geom_text(aes(label= stars), 
            position = position_dodge(0.8), hjust = -1, size=7) +
  
  coord_flip(ylim = c(-20, 155)) +
  theme(axis.text.y = element_text(colour = "black")))

#patchwork
#save as 1100*800
(patch <- (scatter + plot_spacer()) / (p1 + p2 + p3) + 
    plot_layout(heights = c(1,2), guides = "collect") + 
    plot_annotation(tag_levels = "a"))

wrap_elements(panel = patch) +
  labs(tag = "Response (% Change per 50 Mg/ha)") +
  theme(
    plot.tag = element_text(size = 15),
    plot.tag.position = "bottom")
