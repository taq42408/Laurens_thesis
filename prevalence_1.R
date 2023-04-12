flower_meta=read_csv("C:/Users/laure/Desktop/Laurens_thesis/Flower_metadata.csv")
flower_qpcr=read_csv("C:/Users/laure/Desktop/Laurens_thesis/Flower_qPCR_master.csv")%>%
  filter(Content=="Unkn" & !Sample_ID=="NA")
varroa <- read_csv("C:/Users/laure/Desktop/Laurens_thesis/varroa-loads.csv") %>% 
  group_by(Apiary_ID) %>% 
  summarize(varroa_avg = mean(Varroa_Load))

library(dplyr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)

flower_meta_comb=flower_qpcr%>%
  left_join(flower_meta,by=('Sample_ID'))%>%
  mutate(dwv=ifelse(c(SQ>0 & `Melt Temp`<=83 & `Melt Temp`>=81),1,0))%>%
  relocate(dwv,.after=`Melt Temp`)

apiary_size = flower_meta_comb %>% 
  filter(!`# of Colonies`=="NA") %>% 
  group_by(Apiary_ID) %>% 
  summarize(colony_number = mean(`# of Colonies`))
View(apiary_size)

prevalence = table(flower_meta_comb$Apiary_ID, flower_meta_comb$dwv)
prevalence_table = as.data.frame(prop.table(prevalence, margin=1))%>%
  filter(Var2==1)%>%
  rename(Apiary_ID=Var1, dwv_presence=Var2)

natural=read.csv("C:/Users/laure/Desktop/Laurens_thesis/aggregate_values.csv")%>%
  filter(Prelim_Agg_Land_Class=='Natural')%>%
  rename(natural_percent=percent_land_cover)

prevalence_table_2=prevalence_table%>%
  left_join(varroa) %>% 
  left_join(apiary_size) %>% 
  left_join(agriculture)%>%
  left_join(natural)%>%
  select(-X)

v=lm(prevalence_table_2$Freq~prevalence_table_2$varroa_avg)
c=lm(prevalence_table_2$Freq~prevalence_table_2$colony_number)
n=lm(prevalence_table_2$Freq~prevalence_table_2$natural_percent)
a=lm(prevalence_table_2$Freq~prevalence_table_2$agriculture_percent)

plot(prevalence_table_2$Freq~prevalence_table_2$varroa_avg)
abline(v,col="blue")
plot(prevalence_table_2$Freq~prevalence_table_2$colony_number)
abline(c,col="blue")
plot(prevalence_table_2$Freq~prevalence_table_2$natural_percent)
abline(n,col="blue")
plot(prevalence_table_2$Freq~prevalence_table_2$agriculture_percent)
abline(a,col="blue")

agriculture=read.csv("C:/Users/laure/Desktop/Laurens_thesis/aggregrate_values.csv")%>%
  filter(Prelim_Agg_Land_Class=='Agriculture')%>%
  rename(agriculture_percent=percent_land_cover)

# STATS HERE ----
prev = flower_meta_comb %>% 
  group_by(Apiary_ID) %>% 
  summarize(DWV_positives = sum(dwv), total_flowers=n()) %>% 
  mutate(DWV_prev=DWV_positives/total_flowers) %>% 
  mutate(DWV_negatives=(total_flowers-DWV_positives)) %>% 
  left_join(varroa) %>% 
  left_join(apiary_size)
View(prev)

hist(prev$DWV_prev)

model <-glm(cbind(DWV_positives, DWV_negatives)~varroa_avg + colony_number, family=binomial, data=prev)
summary(model)
#successes first, failures second in binomial model 
#As a two-column integer matrix: the first column gives the number of successes and the second the number of failure
#Now we need to check if the model is overdispersed

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(model)
