library(dplyr)
library(tidyverse)
library(ggplot2)
library(forcats)
library(lme4)
library(lmerTest)
library(emmeans)

flower_qpcr <- read_csv("/Users/krd59/Documents/Github/Laurens_thesis/Flower_qPCR_master.csv") %>% 
  filter(Content=="Unkn" & !Sample_ID=="NA")
  
flower_meta <-read_csv("/Users/krd59/Documents/Github/Laurens_thesis/Flower_metadata.csv")%>% 
  rename(Apiary_ID = `Apiary ID`)
 # select(-c("...25","...26","...27", "...28","...29","...30"))
varroa <- read_csv("/Users/krd59/Documents/Github/Laurens_thesis/varroa-loads.csv") %>% 
  group_by(Apiary_ID) %>% 
  summarize(varroa_avg = mean(Varroa_Load))
View(varroa)

treatments <- read_csv("/Users/krd59/Documents/Github/Laurens_thesis/TechTeam-data_2020-ONLY.csv") %>% 
  filter(!`Apiary Treatments per year`=="N/A")
treatments$`Apiary Treatments per year`=as.numeric(treatments$`Apiary Treatments per year`)
View(treatments)

natural=read.csv("/Users/krd59/Documents/Github/Laurens_thesis/aggregate_values.csv")%>%
  filter(Prelim_Agg_Land_Class=='Natural')%>%
  rename(percent_natural=percent_land_cover)


beek <-treatments %>% 
  filter(!Apiary_ID=="NA") %>% 
  group_by(Apiary_ID) %>% 
  summarize(number_treatments = mean(`Apiary Treatments per year`))
unique(beek$Apiary_ID)
  
View(beek)
as.data.frame(flower_qpcr)

flower_meta_comb <- flower_qpcr %>% 
  left_join(flower_meta, by="Sample_ID") %>% 
  mutate(dwv=ifelse(c(SQ>0 & `Melt Temp`<= 83 & `Melt Temp` >= 81), 1, 0)) %>% 
  relocate(dwv, .after=`Melt Temp`) %>% 
  rename(Apiary_ID=Apiary_ID.x) %>% 
  select(Sample_ID, Apiary_ID, Cq, SQ, `Melt Temp`, dwv, Date,`# of Colonies`, `Quadrat #`, `Distance to colonies (m)`, `# of HB Visits`) %>% 
  left_join(varroa)
flower_meta_comb$`Distance to colonies (m)` = as.numeric(flower_meta_comb$`Distance to colonies (m)`)
View(flower_meta_comb)
str(flower_meta_comb)

apiary_size = flower_meta %>% 
  filter(!`# of Colonies`=="NA") %>% 
  group_by(Apiary_ID) %>% 
  summarize(colony_number = mean(`# of Colonies`))
View(apiary_size)
visitation = flower_meta %>% 
  filter(!`# of HB Visits`=="NA") %>% 
  group_by(Apiary_ID) %>% 
  summarize(average_visitation = mean(`# of HB Visits`))
View(apiary_size)

prev = flower_meta_comb %>% 
  group_by(Apiary_ID) %>% 
  summarize(DWV_positives = sum(dwv), total_flowers=n()) %>% 
  mutate(DWV_prev=DWV_positives/total_flowers) %>% 
  mutate(DWV_negatives=(total_flowers-DWV_positives)) %>% 
  left_join(varroa) %>% 
  left_join(apiary_size) %>% 
  left_join(visitation) %>% 
  left_join(beek) %>% 
  left_join(natural)
View(prev)
str(prev)

vc = lm(prev$varroa_avg~prev$colony_number)
plot(prev$DWV_prev~prev$varroa_avg)
plot(prev$DWV_prev~prev$colony_number)
plot(prev$varroa_avg~prev$colony_number)
abline(vc, col="darkred")

plot(flower_meta_comb$dwv~flower_meta_comb$`Distance to colonies (m)`)
plot(flower_meta_comb$dwv~flower_meta_comb$varroa_avg)

prev2 <- prev %>% 
  filter(!number_treatments=="NA")
plot(prev2$varroa_avg~prev2$number_treatments)

round(cor(prev2[, c("varroa_avg","number_treatments", "colony_number", "average_visitation", "percent_natural")]), 3)

model <-glmer(cbind(DWV_positives, DWV_negatives)~varroa_avg+average_visitation+percent_natural+number_treatments+(1|Apiary_ID), family=binomial, data=prev)
summary(model)

# STATS HERE ----
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

#check overdispersion parameter - our model IS significantly overdispersed (p-value is significant);
#to correct for this, we can use an observation-level random effect; each observation has its own p (that is normally distributed)
#in binomials, our parameters are n (total sample size) & p (proportion )
#overdispersion is the presence of greater variability (statistical dispersion) in a data set than would be expected based on a given statistical model.

model <-glmer(cbind(DWV_positives, DWV_negatives)~varroa_avg+colony_number+(1|Apiary_ID), family=binomial, data=prev)
#+Num_colonies as fixed effect
summary(model)
overdisp_fun(model)
#random intercept model 
plot(prev$DWV_prev~prev$average_visitation)
plot(prev$average_visitation~prev$colony_number)

whole.model <-glmer(cbind(DWV_positives, DWV_negatives)~varroa_avg+average_visitation+(1|Apiary_ID), family=binomial, data=prev)
summary(whole.model)

# three reasons for statistics: 1) hypothesis testing; 2) data exploration; 3) prediction modeling
# data exploration - not always compatible for hypothesis testing; you are more likely to select variables that have lower p-values than their average
# so cannot do hypothesis testing on same dataset bc model selection - split dataset for exploratory analysis and half hypothesis
# interactions??? 
# regression to the mean - analogy of choosing 10 fastest runners; they run again and will more likely run slower (the first time they were probably running faster than their personal "true" average)
# check that variables are not correlated with each other: run a correlation matrix
# to deal with co-linearity, we often choose before running the models 
# just because something has a p-value of >0.05, it might still explain some background variation and is fine to keep in the model
# Best statistical practice: for normally distributed data, the min recommendation is 20 points per variable (check on this!)
# run full model, no more than 4 variables 
# biologically relevant interactions - check out 
# include a few potential interactions that you expect to be important, if p-value is insignificant, then you can drop the interaction term 
# in methods, we looked at interactions bw xyz variables, not significant so we dropped these terms


emmeans(model, ~varroa_avg,at=list(varroa_avg= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)), type="response")
model.emp <- emmip(model, ~varroa_avg, at=list(varroa_avg= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)), type="response", CIs = TRUE, plotit = FALSE)
#this at the logit of the probability of the model
#plot data on top of our emmip

# emmeans(model, ~varroa_avg,at=list(varroa_avg= c(1,5,10,15,20)), type="response")
# model.emp <- emmip(model, ~varroa_avg, at=list(varroa_avg= c(1,5,10,15,20)), type="response", CIs = TRUE, plotit = FALSE)
# fewer steps, the more chunky the line looks 

# what is a nested model - every variable present in one model is also in the other model; AIC allows you test bw non-nested values while ANOVA only allows you to compare nested models

ggplot(data=model.emp, aes(x=varroa_avg, y=yvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), alpha=0.1, show.legend = FALSE) + 
  geom_point(data=prev, aes(x=varroa_avg, y=DWV_prev), size=3)+ 
  theme_light() +
  xlab("Average Varroa Load per 100 honey bees") +
  ylab("DWV prevalence on goldenrod flowers") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust=1))

##this model has accounted for the number of colonies! <--- I took it out for simplicity Nov 6 2022, add back in later

emmeans(model, ~colony_number,at=list(colony_number= c(1,5,10,15,20,25,30,35,40,45,50)), type="response")
model.emp.col <- emmip(model, ~colony_number, at=list(colony_number= c(1,5,10,15,20,25,30,35,40,45,50)), type="response", CIs = TRUE, plotit = FALSE)
#this at the logit of the probability of the model
#plot data on top of our emmip

ggplot(data=model.emp.col, aes(x=colony_number, y=yvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), alpha=0.1, show.legend = FALSE) + 
  geom_point(data=prev, aes(x=colony_number, y=DWV_prev), size=3)+ 
  theme_light() +
  xlab("Number of colonies") +
  ylab("DWV prevalence on goldenrod flowers") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust=1))

###### Flower ~ Distance --------
flower_meta_comb <- flower_qpcr %>% 
  left_join(flower_meta, by="Sample_ID") %>% 
  mutate(dwv=ifelse(c(SQ>0 & `Melt Temp`<= 83 & `Melt Temp` >= 81), 1, 0)) %>% 
  relocate(dwv, .after=`Melt Temp`) %>% 
  rename(Apiary_ID=Apiary_ID.x) %>% 
  select(Sample_ID, Apiary_ID, Cq, SQ, `Melt Temp`, dwv, Date,`# of Colonies`, `Quadrat #`, `Distance to colonies (m)`, `# of HB Visits`) %>% 
  left_join(varroa)
flower_meta_comb$`Distance to colonies (m)` = as.numeric(flower_meta_comb$`Distance to colonies (m)`)
View(flower_meta_comb)
str(flower_meta_comb)

model <-glmer(dwv~`Distance to colonies (m)`+(1|Apiary_ID), family=binomial, data=flower_meta_comb)
summary(model)
overdisp_fun(model)
#random intercept model 


emmeans(model, ~`Distance to colonies (m)`, at=list(`Distance to colonies (m)`= c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550)), type="response")
model.emp <- emmip(model, ~`Distance to colonies (m)`, at=list(`Distance to colonies (m)`= c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550)), type="response", CIs = TRUE, plotit = FALSE)
#this at the logit of the probability of the model
#plot data on top of our emmip

#emmeans(model, ~varroa_avg,at=list(varroa_avg= c(1,5,10,15)), type="response")
#model.emp <- emmip(model, ~varroa_avg, at=list(varroa_avg= c(1,5,10,15)), type="response", CIs = TRUE, plotit = FALSE)
#fewer steps, the more chunky the line looks 

ggplot(data=model.emp, aes(x=`Distance to colonies (m)`, y=yvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), alpha=0.1, show.legend = FALSE) + 
  geom_point(data=flower_meta_comb, aes(x=`Distance to colonies (m)`, y=dwv), size=3)+ 
  theme_light() +
  xlab("Distance from nearest honey bee colony") +
  ylab("DWV prevalence on goldenrod flowers") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust=1))
  