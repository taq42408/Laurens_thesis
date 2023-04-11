flower_meta=read_csv("C:/Users/laure/Desktop/Laurens_thesis/Flower_metadata.csv")%>%
  select(-c("...25","...26","...27","...28","...29","...30"))
flower_qpcr=read_csv("C:/Users/laure/Desktop/Laurens_thesis/Flower_qPCR_master.csv")
library(dplyr)
library(tidyverse)
library(ggplot2)
