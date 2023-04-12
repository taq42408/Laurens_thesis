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
s=lm(prevalence_table_2$varroa_avg~prevalence_table_2$colony_number)

plot(prevalence_table_2$varroa_avg~prevalence_table_2$colony_number)
abline(s,col="blue")
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

#check overdispersion parameter - our model IS significantly overdispersed (p-value is significant);
#to correct for this, we can use an observation-level random effect; each observation has its own p (that is normally distributed)
#in binomials, our parameters are n (total sample size) & p (proportion )
#overdispersion is the presence of greater variability (statistical dispersion) in a data set than would be expected based on a given statistical model.
#in binomial/poisson, the mean and variance are related to each other (binomial - probabilityy of each sample turning up positive is exactly the same)
#So for 5 sites, we have all the same treatments, we assume the probability of positive is the same, and mean and variance are assumed to be the same 
# but there are inter-site differences that are not the treatments that influence the probability, so we need to account for these site to site differences w a random effect of site
# random effect removes the link bw mean and variance 

model <-glmer(cbind(DWV_positives, DWV_negatives)~varroa_avg+colony_number+(1|Apiary_ID), family=binomial, data=prev)
#+Num_colonies as fixed effect
summary(model)
overdisp_fun(model)
#random intercept model 

emmeans(model, ~varroa_avg,at=list(varroa_avg= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)), type="response")
model.emp <- emmip(model, ~varroa_avg, at=list(varroa_avg= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)), type="response", CIs = TRUE, plotit = FALSE)
#this at the logit of the probability of the model
#plot data on top of our emmip

#emmeans(model, ~varroa_avg,at=list(varroa_avg= c(1,5,10,15)), type="response")
#model.emp <- emmip(model, ~varroa_avg, at=list(varroa_avg= c(1,5,10,15)), type="response", CIs = TRUE, plotit = FALSE)
#fewer steps, the more chunky the line looks 

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
