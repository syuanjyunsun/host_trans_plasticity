#### R scripts for Transgenerational plasticity in a zooplankton in response to temperature elevation and parasitism ####

data=read.csv("data_for_host_trans_plasticity.csv",header=T)

#a subset of hosts that survived during processing
data2=data[data$earlydeath=="0",]

library(lme4)
library(car)
library(emmeans)

#model selection based on AICc
# load the MUMIn package 
require(MuMIn)

####host traits for initial infection process####
model=glmer(log(Body.size+1)~F1temp+(1|source),gaussian,data=data)
Anova(model,type=3)

#normality testing for gut.resistance
datagut.resistance=data[data$F1parasite=="infection",]
shapiro.test(datagut.resistance$gut.resistance)

model1.1=glmer(gut.resistance~(F0temp)*(F0parasite)*F1temp+mean.foregut.width+(1|source),gaussian,data=data[data$F1parasite=="infection",])
model1.2=glmer(gut.resistance~(F0temp+F0parasite)*F1temp+F0temp*F0parasite+mean.foregut.width+(1|source),gaussian,data=data[data$F1parasite=="infection",])
model1.3=glmer(gut.resistance~(F0temp)*F1temp+F0temp*F0parasite+mean.foregut.width+(1|source),gaussian,data=data[data$F1parasite=="infection",])
model1.4=glmer(gut.resistance~(F0temp)+F1temp+F0temp*F0parasite+mean.foregut.width+(1|source),gaussian,data=data[data$F1parasite=="infection",])
model1.5=glmer(gut.resistance~(F0temp)+(F0parasite)+F1temp+mean.foregut.width+(1|source),gaussian,data=data[data$F1parasite=="infection",])
out.put<-model.sel(model1.1,model1.2,model1.3,model1.4,model1.5)
out.put

Anova(model1.5,type=3)

#normality testing for Hemocytes.per.spore
dataHemocytes.per.spore=data[data$F1parasite=="infection",]
shapiro.test(dataHemocytes.per.spore$Hemocytes.per.spore)

####parasite traits for initial infection process####
model2.1=glmer(log(Hemocytes.per.spore+1)~F0temp*F1temp*F0parasite+(1|source),gaussian,data=data[data$F1parasite=="infection",])
model2.2=glmer(log(Hemocytes.per.spore+1)~F0temp*F1temp+(F0temp+F1temp)*F0parasite+(1|source),gaussian,data=data[data$F1parasite=="infection",])
model2.3=glmer(log(Hemocytes.per.spore+1)~(F0temp+F1temp)*F0parasite+(1|source),gaussian,data=data[data$F1parasite=="infection",])
model2.4=glmer(log(Hemocytes.per.spore+1)~F0temp+(F1temp)*F0parasite+(1|source),gaussian,data=data[data$F1parasite=="infection",])
model2.5=glmer(log(Hemocytes.per.spore+1)~F0temp+(F1temp)+F0parasite+(1|source),gaussian,data=data[data$F1parasite=="infection",])
out.put<-model.sel(model2.1,model2.2,model2.3,model2.4,model2.5)
out.put

Anova(model2.5,type=3)

####age at first reproduction####
model3.1=glmer(ageatfirst~F0temp*F1temp*F0parasite*F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.2=glmer(ageatfirst~(F0temp+F1temp)*F0parasite*F1parasite+F0temp*F1temp*(F0parasite+F1parasite)+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.3=glmer(ageatfirst~(F0temp)*F0parasite*F1parasite+F0temp*F1temp*(F0parasite+F1parasite)+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.4=glmer(ageatfirst~(F0temp)*F0parasite*F1parasite+F0temp*F1temp*(F1parasite)+F0parasite*F1temp+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.5=glmer(ageatfirst~(F0temp)*F0parasite*F1parasite+F0temp*(F1temp+F1parasite)+F1temp*F1parasite+F0parasite*F1temp+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.6=glmer(ageatfirst~F0temp*F0parasite+(F0temp+F0parasite)*F1parasite+F0temp*(F1temp+F1parasite)+F1temp*F1parasite+F0parasite*F1temp+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.7=glmer(ageatfirst~F0temp*F0parasite+(F0temp+F0parasite)*F1parasite+F0temp*(F1temp+F1parasite)+F1temp+F1parasite+F0parasite*F1temp+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.8=glmer(ageatfirst~F0temp*F0parasite+(F0temp)*F1parasite+F0temp*(F1temp+F1parasite)+F1temp+F1parasite+F0parasite*F1temp+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.9=glmer(ageatfirst~F0temp*F0parasite+F0temp*(F1temp)+F1temp+F1parasite+F0parasite*F1temp+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.10=glmer(ageatfirst~F0temp*F0parasite+F0temp*(F1temp)+F1temp+F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.11=glmer(ageatfirst~F0temp*F0parasite+F1temp+F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model3.12=glmer(ageatfirst~F0temp+F1temp+F0parasite+F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
out.put<-model.sel(model3.1,model3.2,model3.3,model3.4,model3.5,model3.6,model3.7,model3.8,model3.9,model3.10,model3.11,model3.12)
out.put

Anova(model3.12,type=3)

####broodsize of first reproduction####
model4.1=glmer(firstclutch~F0temp*F1temp*F0parasite*F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.2=glmer(firstclutch~F0temp*F1temp*(F0parasite+F1parasite)+(F0temp+F1temp)*F0parasite*F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.3=glmer(firstclutch~F0temp*F1temp*(F0parasite)+(F0temp+F1temp)*F0parasite*F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.4=glmer(firstclutch~F0temp*(F1temp+F0parasite)+(F0temp+F1temp)*F0parasite*F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.5=glmer(firstclutch~F0temp*(F1temp+F0parasite)+(F0temp)*F0parasite*F1parasite+F1temp*(F0parasite+F1parasite)+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.6=glmer(firstclutch~F0temp*(F1temp+F0parasite)+(F0temp)*(F0parasite+F1parasite)+F1temp*(F0parasite+F1parasite)+F0parasite*F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.7=glmer(firstclutch~F0temp*(F1temp)+(F0temp)*(F1parasite)+F1temp*(F0parasite+F1parasite)+F0parasite*F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.8=glmer(firstclutch~F0temp*(F1temp)+F1temp*(F0parasite+F1parasite)+F0parasite*F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.9=glmer(firstclutch~F0temp*(F1temp)+F1temp*(F1parasite)+F0parasite*F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.10=glmer(firstclutch~F0temp*(F1temp)+F1temp*(F1parasite)+F0parasite+F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.11=glmer(firstclutch~F0temp*(F1temp)+F1temp+(F1parasite)+F0parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
model4.12=glmer(firstclutch~F0temp+F1temp+F0parasite+F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
out.put<-model.sel(model4.1,model4.2,model4.3,model4.4,model4.5,model4.6,model4.7,model4.8,model4.9,model4.10,model4.11,model4.12)
out.put

anova(model4.11,model4.12)
#model4.11 and model4.12 are equally viable, but non-significant interaction was preferably removed from model4.11, making model4.12 the final model. 

Anova(model4.12,type=3)

model5.1=coxme(Surv(lifespan)~F0temp*F1temp*F0parasite*F1parasite+(1|source),data=data2[data2$infection.status!="exposed-uninfected",])

#Significant four-way interactions were detected.
Anova(model5.1,type=3)

#evaluating the proportional hazards assumption
library("survival")
library("survminer")
km <- survfit(Surv(lifespan)~F1parasiteF1temp,data=data2[data2$infection.status!="exposed-uninfected",])
ggsurvplot(km, fun = "cloglog")
#the plotting of log-log Kaplan Meier survival estimates against time between treatments are reasonably parallel.

#the Schoenfeld residuals examining model fit also yields flat values centered around zero.
zph <- cox.zph(model5.1)
par(mfrow = c(1, 2))
plot(zph, var = 1)
plot(zph, var = 2)

#the goodness-of-fit test also suggests a correlation of zero between the Schoenfeld residuals and survival time, indicating the model met the assumption. 
cox.zph(model5.1)

#post-hoc comparisons
posthoc=emmeans(model5.1,  ~  F1temp+F1parasite|F0temp+F0parasite, adjust="tukey")
pairs(posthoc)

####lifetime fecundity####
model6.1=glmer.nb(totalbaby~F0temp*F1temp*F0parasite*F1parasite+(1|source),data=data2[data2$infection.status!="exposed-uninfected",])

#Significant four-way interactions were detected.
Anova(model6.1,type=3)

#post-hoc comparisons
posthoc=emmeans(model6.1,  ~  F1temp+F1parasite|F0temp+F0parasite, adjust="tukey")
pairs(posthoc)

#considering immune responses in predicting fecundity
model6.1.1=glmer.nb(totalbaby~F0parasite*(F1temp*F0temp)*(Hemocytes.per.spore+gut.resistance)+(1|source),data=data2[data2$infection=="1",])
model6.1.2=glmer.nb(totalbaby~F0parasite*(F1temp*F0temp)*(Hemocytes.per.spore)+F0parasite*(F1temp+F0temp)*gut.resistance+F1temp*F0temp*(F0parasite+gut.resistance)+(1|source),data=data2[data2$infection=="1",])
model6.1.3=glmer.nb(totalbaby~F0parasite*(F1temp+F0temp)*(Hemocytes.per.spore)+(Hemocytes.per.spore+F0parasite)*(F1temp*F0temp)+F0parasite*(F1temp+F0temp)*gut.resistance+F1temp*F0temp*(F0parasite+gut.resistance)+(1|source),data=data2[data2$infection=="1",])
model6.1.4=glmer.nb(totalbaby~F0parasite*(F1temp+F0temp)*(Hemocytes.per.spore)+(Hemocytes.per.spore+F0parasite)*(F1temp*F0temp)+F0parasite*(F1temp)*gut.resistance+F1temp*F0temp*(F0parasite+gut.resistance)+(1|source),data=data2[data2$infection=="1",])
model6.1.5=glmer.nb(totalbaby~F0parasite*(F1temp+F0temp)*(Hemocytes.per.spore)+(Hemocytes.per.spore+F0parasite)*(F1temp*F0temp)+F0parasite*(F1temp)*gut.resistance+F1temp*F0temp*(F0parasite)+gut.resistance*F0temp+(1|source),data=data2[data2$infection=="1",])
model6.1.6=glmer.nb(totalbaby~F0parasite*(F1temp+F0temp)*(Hemocytes.per.spore)+(Hemocytes.per.spore+F0parasite)*(F1temp*F0temp)+(F0parasite+F1temp)*gut.resistance+F1temp*F0temp*(F0parasite)+gut.resistance*F0temp+(1|source),data=data2[data2$infection=="1",])
model6.1.7=glmer.nb(totalbaby~F0parasite*(F1temp+F0temp)*(Hemocytes.per.spore)+(Hemocytes.per.spore)*(F1temp*F0temp)+(F0parasite+F1temp)*gut.resistance+F1temp*(F0temp+F0parasite)+gut.resistance*F0temp+(1|source),data=data2[data2$infection=="1",])
model6.1.8=glmer.nb(totalbaby~F0parasite*(F1temp+F0temp)*(Hemocytes.per.spore)+(F0parasite+F1temp)*gut.resistance+F1temp*(F0temp+F0parasite)+gut.resistance*F0temp+(1|source),data=data2[data2$infection=="1",])
model6.1.9=glmer.nb(totalbaby~F0parasite*(F1temp+F0temp)*(Hemocytes.per.spore)+(F1temp)*gut.resistance+F1temp*(F0temp+F0parasite)+gut.resistance*F0temp+(1|source),data=data2[data2$infection=="1",])
model6.1.10=glmer.nb(totalbaby~F0parasite*(F1temp+F0temp)*(Hemocytes.per.spore)+(F1temp)+gut.resistance+F1temp*(F0temp+F0parasite)+gut.resistance*F0temp+(1|source),data=data2[data2$infection=="1",])
model6.1.11=glmer.nb(totalbaby~F0parasite*(F1temp+F0temp)*(Hemocytes.per.spore)+(F1temp)+gut.resistance+F1temp*(F0temp+F0parasite)+gut.resistance+F0temp+(1|source),data=data2[data2$infection=="1",])
out.put<-model.sel(model6.1.1,model6.1.2,model6.1.3,model6.1.4,model6.1.5,model6.1.6,model6.1.7,model6.1.8,model6.1.9,model6.1.10,model6.1.11)
out.put

Anova(model6.1.11,type=3)

####probability of infection####
model7.1=glmer(infection~F0temp*F1temp*F0parasite+(1|source),binomial,data=data2)
model7.2=glmer(infection~(F0temp+F1temp)*F0parasite+F0temp*F1temp+(1|source),binomial,data=data2)
model7.3=glmer(infection~(F0temp)*F0parasite+F0temp*F1temp+(1|source),binomial,data=data2)
model7.4=glmer(infection~(F0temp)*F0parasite+F0temp+F1temp+(1|source),binomial,data=data2)
model7.5=glmer(infection~(F0temp)+F0parasite+F0temp+F1temp+(1|source),binomial,data=data2)
out.put<-model.sel(model7.1,model7.2,model7.3,model7.4,model7.5)
out.put

Anova(model7.5,type=3)

#normality testing for Hemocytes.per.spore
datamaturetrans=data2[data2$infection=="1",]
shapiro.test(datamaturetrans$maturetrans)

####spore production####
model8.1=glmer(log(maturetrans+1)~F0temp*F1temp*F0parasite+gut.resistance+Hemocytes.per.spore+(1|source),gaussian,data=data2[data2$infection=="1",])
model8.2=glmer(log(maturetrans+1)~F0temp*F1temp+(F0temp+F1temp)*F0parasite+gut.resistance+Hemocytes.per.spore+(1|source),gaussian,data=data2[data2$infection=="1",])
model8.3=glmer(log(maturetrans+1)~F0temp+F1temp+(F0temp+F1temp)*F0parasite+gut.resistance+Hemocytes.per.spore+(1|source),gaussian,data=data2[data2$infection=="1",])
model8.4=glmer(log(maturetrans+1)~F0temp+F1temp+(F1temp)*F0parasite+gut.resistance+Hemocytes.per.spore+(1|source),gaussian,data=data2[data2$infection=="1",])
model8.5=glmer(log(maturetrans+1)~F0temp+F1temp+(F1temp)+F0parasite+gut.resistance+Hemocytes.per.spore+(1|source),gaussian,data=data2[data2$infection=="1",])
out.put<-model.sel(model8.1,model8.2,model8.3,model8.4,model8.5)
out.put

Anova(model8.5,type=3)


