library(dplyr)
library(tidyverse)
library(ggplot2)
library(forcats)
library(lme4)
library(lmerTest)
library(emmeans)

flower_qpcr <- read_csv("/Users/krd59/Documents/Github/Laurens_thesis/Flower_qPCR_master.csv") %>% 
  filter(Content=="Unkn" & !Sample_ID=="NA")
  
flower_meta <-read_csv("/Users/krd59/Documents/Github/Laurens_thesis/Flower_metadata.csv")
#%>% 
 # select(-c("...25","...26","...27", "...28","...29","...30"))
varroa <- read_csv("/Users/krd59/Documents/Github/Laurens_thesis/varroa-loads.csv") %>% 
  group_by(Apiary_ID) %>% 
  summarize(varroa_avg = mean(Varroa_Load))
  
View(flower_qpcr)
as.data.frame(flower_qpcr)

flower_meta_comb <- flower_qpcr %>% 
  left_join(flower_meta, by="Sample_ID") %>% 
  mutate(dwv=ifelse(c(SQ>0 & `Melt Temp`<= 83 & `Melt Temp` >= 81), 1, 0)) %>% 
  relocate(dwv, .after=`Melt Temp`)
View(flower_meta_comb)

prev = flower_meta_comb %>% 
  group_by(Apiary_ID) %>% 
  summarize(DWV_positives = sum(dwv), total_flowers=n()) %>% 
  mutate(DWV_prev=DWV_positives/total_flowers) %>% 
  mutate(DWV_negatives=(total_flowers-DWV_positives)) %>% 
  left_join(varroa) %>% 
  left_join(apiary_size)

apiary_size = flower_meta_comb %>% 
  filter(!`# of Colonies`=="NA") %>% 
  group_by(Apiary_ID) %>% 
  summarize(colony_number = mean(`# of Colonies`))
View(apiary_size)

plot(prev$DWV_prev~prev$varroa_avg)
plot(prev$DWV_prev~prev$colony_number)
  