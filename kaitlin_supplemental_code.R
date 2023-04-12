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
  
View(flower_qpcr)
as.data.frame(flower_qpcr)

flower_meta_comb <- flower_qpcr %>% 
  left_join(flower_meta, by="Sample_ID") %>% 
  mutate(dwv=ifelse(c(SQ>0 & `Melt Temp`<= 83 & `Melt Temp` >= 81), 1, 0)) %>% 
  relocate(dwv, .after=`Melt Temp`) %>% 
  rename(Apiary_ID=Apiary_ID.x) %>% 
  select(Sample_ID, Apiary_ID, Cq, SQ, `Melt Temp`, dwv, Date,`# of Colonies`, `Quadrat #`, `Distance to colonies (m)`, `# of HB Visits`) %>% 
  left_join(varroa)

View(flower_meta_comb)
str(flower_meta_comb)
flower_meta_comb$`Distance to colonies (m)` = as.numeric(flower_meta_comb$`Distance to colonies (m)`)

apiary_size = flower_meta %>% 
  filter(!`# of Colonies`=="NA") %>% 
  group_by(Apiary_ID) %>% 
  summarize(colony_number = mean(`# of Colonies`))
View(apiary_size)

prev = flower_meta_comb %>% 
  group_by(Apiary_ID) %>% 
  summarize(DWV_positives = sum(dwv), total_flowers=n()) %>% 
  mutate(DWV_prev=DWV_positives/total_flowers) %>% 
  mutate(DWV_negatives=(total_flowers-DWV_positives)) %>% 
  left_join(varroa) %>% 
  left_join(apiary_size)
View(prev)

vc = lm(prev$varroa_avg~prev$colony_number)
plot(prev$DWV_prev~prev$varroa_avg)
plot(prev$DWV_prev~prev$colony_number)
plot(prev$varroa_avg~prev$colony_number)
abline(vc, col="darkred")

plot(flower_meta_comb$dwv~flower_meta_comb$`Distance to colonies (m)`)
plot(flower_meta_comb$dwv~flower_meta_comb$varroa_avg)


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



  