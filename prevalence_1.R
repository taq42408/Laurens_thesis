library(dplyr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)

flower_meta=read_csv("C:/Users/laure/Desktop/Laurens_thesis/Flower_metadata.csv")
flower_qpcr=read_csv("C:/Users/laure/Desktop/Laurens_thesis/Flower_qPCR_master.csv")%>%
  filter(Content=="Unkn" & !Sample_ID=="NA")
varroa <- read_csv("C:/Users/laure/Desktop/Laurens_thesis/varroa-loads.csv") %>% 
  group_by(Apiary_ID) %>% 
  summarize(varroa_avg = mean(Varroa_Load))


flower_meta_comb=flower_qpcr%>%
  left_join(flower_meta,by=('Sample_ID'))%>%
  mutate(dwv=ifelse(c(SQ>0 & `Melt Temp`<=83 & `Melt Temp`>=81),1,0))%>%
  relocate(dwv,.after=`Melt Temp`) %>% 
  rename(Apiary_ID=Apiary_ID.x)

flower_meta_comb$quadrat_date <- paste(flower_meta_comb$`Quadrat #`,
                                       flower_meta_comb$`Date (numeric)`,
                                       sep="_"
                                       )
flower_meta_comb$`# of HB Visits` <- as.numeric(flower_meta_comb$`# of HB Visits`)
flower_meta_comb$`Distance to colonies (m)`<-as.numeric(flower_meta_comb$`Distance to colonies (m)`)

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
  rename(percent_natural=percent_land_cover)

agriculture=read.csv("C:/Users/laure/Desktop/Laurens_thesis/aggregate_values.csv")%>%
  filter(Prelim_Agg_Land_Class=='Agriculture')%>%
  rename(agriculture_percent=percent_land_cover)

prevalence_table_2=prevalence_table%>%
  left_join(varroa) %>% 
  left_join(apiary_size) %>% 
  left_join(natural)%>%
  select(-X)

v=lm(prevalence_table_2$Freq~prevalence_table_2$varroa_avg)
c=lm(prevalence_table_2$Freq~prevalence_table_2$colony_number)
n=lm(prevalence_table_2$Freq~prevalence_table_2$percent_natural)
s=lm(prevalence_table_2$varroa_avg~prevalence_table_2$colony_number)

plot(prevalence_table_2$varroa_avg~prevalence_table_2$colony_number)
abline(s,col="blue")
plot(prevalence_table_2$Freq~prevalence_table_2$varroa_avg)
abline(v,col="blue")
plot(prevalence_table_2$Freq~prevalence_table_2$colony_number)
abline(c,col="blue")
plot(prevalence_table_2$Freq~prevalence_table_2$percent_natural)
abline(n,col="blue")

honeybee_visits = flower_meta_comb %>% 
  filter(!`# of HB Visits`=="NA") %>% 
  group_by(Apiary_ID) %>% 
  summarize(average_visitation = mean(`# of HB Visits`))

treatments <- read_csv("C:/Users/laure/Desktop/Laurens_thesis/TechTeam-data_2020-ONLY.csv") %>% 
  filter(!`Apiary Treatments per year`=="N/A")
treatments$`Apiary Treatments per year`=as.numeric(treatments$`Apiary Treatments per year`)
beek <-treatments %>% 
  filter(!Apiary_ID=="NA") %>% 
  group_by(Apiary_ID) %>% 
  summarize(number_treatments = mean(`Apiary Treatments per year`))
View(treatments)

# STATS HERE ----
prev = flower_meta_comb %>% 
  group_by(Apiary_ID) %>% 
  summarize(DWV_positives = sum(dwv), total_flowers=n()) %>% 
  mutate(DWV_prev=DWV_positives/total_flowers) %>% 
  mutate(DWV_negatives=(total_flowers-DWV_positives)) %>% 
  left_join(varroa) %>% 
  left_join(apiary_size)%>%
  left_join(honeybee_visits)%>%
  left_join(natural)%>%
  left_join(beek) %>% 
  select(-X)
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

#glmer has fixed and random effects, glm just has fixed effects. g means you are not assuming the data is normally distributed 

prev2=prev %>%
  mutate(average_visitation_min=average_visitation/5)

model <-glmer(cbind(DWV_positives, DWV_negatives)~varroa_avg+average_visitation_min+percent_natural+(1|Apiary_ID), family=binomial, data=prev2)
#+Num_colonies as fixed effect
summary(model)
overdisp_fun(model)
#random intercept model 
round(cor(prev2[, c("varroa_avg", "colony_number", "average_visitation_min", "percent_natural")]), 3)

library(emmeans)
emmeans(model, ~varroa_avg,at=list(varroa_avg= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)), type="response")
model.emp <- emmip(model, ~varroa_avg, at=list(varroa_avg= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)), type="response", CIs = TRUE, plotit = FALSE)
#this at the logit of the probability of the model
#plot data on top of our emmip

#emmeans(model, ~varroa_avg,at=list(varroa_avg= c(1,5,10,15)), type="response")
#model.emp <- emmip(model, ~varroa_avg, at=list(varroa_avg= c(1,5,10,15)), type="response", CIs = TRUE, plotit = FALSE)
#fewer steps, the more chunky the line looks 
#emmeans is determining the predictive prevalence of DWV at each site, which we then draw a line through with the model function

ggplot(data=model.emp, aes(x=varroa_avg, y=yvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), alpha=0.1, show.legend = FALSE) + 
  geom_point(data=prev2, aes(x=varroa_avg, y=DWV_prev), size=3)+ 
  theme_light() +
  xlab("Mean Varroa load per apiary") +
  ylab("DWV prevalence on goldenrod flowers") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust=1))

emmeans(model, ~colony_number,at=list(colony_number= c(1,5,10,15,20,25,30,35,40,45,50)), type="response")
model.emp.col <- emmip(model, ~colony_number, at=list(colony_number= c(1,5,10,15,20,25,30,35,40,45,50)), type="response", CIs = TRUE, plotit = FALSE)
#this at the logit of the probability of the model
#plot data on top of our emmip

ggplot(data=model.emp.col, aes(x=colony_number, y=yvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), alpha=0.1, show.legend = FALSE) + 
  geom_point(data=prev2, aes(x=colony_number, y=DWV_prev), size=3)+ 
  theme_light() +
  xlab("Number of colonies") +
  ylab("DWV prevalence on goldenrod flowers") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust=1))

emmeans(model, ~average_visitation_min,at=list(average_visitation_min= c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)), type="response")
model.emp.hb <- emmip(model, ~average_visitation_min, at=list(average_visitation_min= c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)), type="response", CIs = TRUE, plotit = FALSE)

ggplot(data=model.emp.hb, aes(x=average_visitation_min, y=yvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), alpha=0.1, show.legend = FALSE) + 
  geom_point(data=prev2, aes(x=average_visitation_min, y=DWV_prev), size=3)+ 
  theme_light() +
  xlab("Honey bee visits/min to goldenrod flowers") +
  ylab("DWV prevalence on goldenrod flowers") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust=1))

emmeans(model, ~percent_natural,at=list(percent_natural= c(0,10,20,30,40,50,60,70,80,90,100)), type="response")
model.emp.nat <- emmip(model, ~percent_natural, at=list(percent_natural= c(0,10,20,30,40,50,60,70,80,90,100)), type="response", CIs = TRUE, plotit = FALSE)
summary(model.emp.nat)

ggplot(data=model.emp.nat, aes(x=percent_natural, y=yvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), alpha=0.1, show.legend = FALSE) + 
  geom_point(data=prev2, aes(x=percent_natural, y=DWV_prev), size=3)+ 
  theme_light() +
  xlab("Percent Natural Land Cover") +
  ylab("DWV prevalence on goldenrod flowers") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust=1))

###### Flower ~ Distance --------

flower_meta_comb$`Distance to colonies (m)` = as.numeric(flower_meta_comb$`Distance to colonies (m)`)
View(flower_meta_comb)
str(flower_meta_comb)

library(glmmTMB)
dist.model <-glmmTMB(dwv~`Distance to colonies (m)`+(1|Apiary_ID/quadrat_date), family=binomial, data=flower_meta_comb)
summary(dist.model)
overdisp_fun(dist.model)
#this model is not overdispersed - we don't need to add a random effect

#glmer is running the statistics to determine significance
#emmeans is determining the predictive chance of one flower having DWV at each distance, which we then draw a line through with the model function
emmeans(dist.model, ~`Distance to colonies (m)`, at=list(`Distance to colonies (m)`= c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550)), type="response")
model.dist.emp <- emmip(dist.model, ~`Distance to colonies (m)`, at=list(`Distance to colonies (m)`= c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550)), type="response", CIs = TRUE, plotit = FALSE)

ggplot(data=model.dist.emp, aes(x=`Distance to colonies (m)`, y=yvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), alpha=0.1, show.legend = FALSE) + 
  geom_point(data=flower_meta_comb, aes(x=`Distance to colonies (m)`, y=dwv), size=3)+ 
  theme_light() +
  xlab("Distance from apiary (m)") +
  ylab("Probability of DWV presence on goldenrod flowers") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust=1))


#visitation on y axis (group by quadrat_date and apiary), distance on x axis
#try linear model first, look at histogram to see if it's normal (lmer or glmer)
#use old spreadsheet with one line per quadrat, need to add quadrat_date

#honeybee_visits = flower_meta_comb %>% 
#filter(!`# of HB Visits`=="NA") %>% 
  #group_by(Apiary_ID) %>% 
  #summarize(average_visitation = mean(`# of HB Visits`))

#probability of DWV on y, honeybee visitation to quadrat on x

#regress visitation and distance - plot
#look to see if it's normal with hist
# regression statistics, emmeans
#if they're not correlated, include visitation as fixed effect on 195

#Correlation between visitation and distance - no correlation

plot(visitation.distance.df$quadrat_visitation~visitation.distance.df$quadrat_dist)
hist(visitation.distance.df$quadrat_dist)
hist(visitation.distance.df$quadrat_visitation)
     
visitation.distance.df <- flower_meta_comb %>% 
  group_by(Apiary_ID, quadrat_date) %>% 
  summarize(quadrat_visitation = mean(`# of HB Visits`), quadrat_dist = mean(`Distance to colonies (m)`))%>%
  filter(!quadrat_visitation=="NA", !quadrat_dist=="NA")
View(flower_meta_comb)
View(visitation.distance.df)

dist_visit=glmmTMB(quadrat_visitation~quadrat_dist+(1|Apiary_ID), family=nbinom2, data=visitation.distance.df)
summary(dist_visit)
overdisp_fun(dist_visit)

#This was Wee Hao trying to determine whether apiary variation is driving visitation differences as opposed to variation in distance. It is- so we can't show a correlation between visitation and distance with our data

palette = rainbow(15)
names(palette) <- unique(visitation.distance.df$Apiary_ID)
pch = 1:15
names(pch) <- unique(visitation.distance.df$Apiary_ID)

plot(visitation.distance.df$quadrat_visitation~visitation.distance.df$quadrat_dist, col=palette[visitation.distance.df$Apiary_ID],
     pch=pch[visitation.distance.df$Apiary_ID], cex=1.5)

#HB Visitation DWV Data

hb.model <-glmmTMB(dwv~`# of HB Visits`+(1|Apiary_ID/quadrat_date), family=binomial, data=flower_meta_comb)
summary(hb.model)
overdisp_fun(hb.model)

emmeans(hb.model, ~`# of HB Visits`, at=list(`# of HB Visits`= c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60)), type="response")
model.visit <- emmip(hb.model, ~`# of HB Visits`, at=list(`# of HB Visits`= c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60)), type="response", CIs = TRUE, plotit = FALSE)

ggplot(data=model.visit, aes(x=`# of HB Visits`, y=yvar)) +
  geom_line() +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), alpha=0.1, show.legend = FALSE) + 
  geom_point(data=flower_meta_comb, aes(x=`# of HB Visits`, y=dwv), size=3)+ 
  theme_light() +
  xlab("Number of honeybee visits per 5 minute interval") +
  ylab("Probability of DWV presence on goldenrod flowers") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust=1))

#Descriptive DWV prevalence graph

prev3=prev2 %>%
  mutate(sd=sd(DWV_prev))


ggplot(data=prev3, aes(x=fct_reorder(Apiary_ID, DWV_prev), y= DWV_prev))+
  geom_col()+
  scale_x_discrete(labels=c('Apiary 1', 'Apiary 2', 'Apiary 3', 'Apiary 4','Apiary 5','Apiary 6','Apiary 7','Apiary 8','Apiary 9','Apiary 10','Apiary 11','Apiary 12','Apiary 13','Apiary 14','Apiary 15'))+
  theme_light()+
  xlab("Apiary ID")+
  ylab("DWV prevalence on goldenrod flowers")+
  geom_errorbar(aes(ymin = DWV_prev-sd, ymax = DWV_prev+sd),width=.2)+
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust=1))

#average prevalence
mean(prev2$DWV_prev)
