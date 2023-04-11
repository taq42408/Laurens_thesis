flower_meta=read_csv("C:/Users/laure/Desktop/Laurens_thesis/Flower_metadata.csv")%>%
  select(-c("...25","...26","...27","...28","...29","...30"))
flower_qpcr=read_csv("C:/Users/laure/Desktop/Laurens_thesis/Flower_qPCR_master.csv")%>%
  filter(Content=="Unkn" & !Sample_ID=="NA")
library(dplyr)
library(tidyverse)
library(ggplot2)
flower_meta_comb=flower_qpcr%>%
  left_join(flower_meta,by=('Sample_ID'))%>%
  mutate(dwv=ifelse(c(SQ>0 & `Melt Temp`<=83 & `Melt Temp`>=81),1,0))%>%
  relocate(dwv,.after=`Melt Temp`)
prevalence = table(flower_meta_comb$Apiary_ID, flower_meta_comb$dwv)
prevalence_table = prop.table(prevalence, margin=1)
prevalence_table
