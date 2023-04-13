
library(tidyverse)
library(dplyr)
library(ggplot2)

Land_Cover = read.csv("buffer_landscape.csv")
crops_codes = read.csv("crop_codes.csv")
Land_Cover2 = Land_Cover %>%
  select(-c("fid","Beekeeper.ID","Operation.Scale", "Apiary.1.Address", "Apiary.1.City", "Apiary.1.State", 
            "Apiary.1.ZIP", "Apiary.1.County", "Apiary.Address", "Before.Tech.Sampling", 
            "After.Tech.Sampling"))
View(Land_Cover2)

## We want to clean up this code - calculate total area in the buffer 
# around each site (RowSums), and reclassify land cover types from "CDL_1" to "corn"

Land_Cover_Cleaned = Land_Cover2 %>%
  mutate(total_area = rowSums(Land_Cover2[,4:44])) %>%
  mutate(across(CDL_1:CDL_242)/total_area)
  
Land_Cover_Cleaned_longer = pivot_longer(Land_Cover_Cleaned, cols = CDL_1:CDL_242, 
                                         names_to = 'Crop_Code', values_to = 'Value') %>%
  mutate(Value = Value*100) %>%
  left_join(crops_codes, by="Crop_Code")

aggregate_values = 
  filter(Land_Cover_Cleaned_longer,Prelim_Agg_Land_Class %in% 
           c('Agriculture','Developed','Natural')) %>%
  group_by(Apiary.1.Name, Prelim_Agg_Land_Class) %>%
  summarize(percent_land_cover = sum(Value)) 

ggplot(data = aggregate_values)+
  geom_histogram(mapping = aes(x = percent_land_cover, 
                               fill=Prelim_Agg_Land_Class))+
  facet_wrap(~Prelim_Agg_Land_Class)+
  theme(legend.position="none")+
  ylab("Number of Sites")+
  xlab("Percent Land Cover")+
  ggtitle("Distribution of Percent Land Cover by Preliminary Aggregate Land Class")


# Made a ggplot
#aggregate_values %>%
#ggplot(aes(x=Apiary.1.Name, y=percent_land_cover, fill=Prelim_Agg_Land_Class))+
  #geom_bar(stat="identity", width = 0.5, position="dodge")+
  #theme_light()+
  #theme(axis.text.x = element_text(angle=75, size=12, hjust=1))+
  #facet_wrap(~Prelim_Agg_Land_Class)

write.csv(aggregate_values, file="aggregate_values.csv")
