#### R scripts for Transgenerational plasticity in a zooplankton in response to temperature elevation and parasitism ####

data=read.csv("data_for_host_trans_plasticity.csv",header=T)

#a subset of hosts that survived during processing
data2=data[data$earlydeath=="0",]

####host traits for initial infection process####
model=glmer(log(Body.size+1)~F1temp+(1|source),gaussian,data=data)
Anova(model,type=3)

model=glmer(gut.resistance~(F0temp)+(F0parasite)+mean.foregut.width+F1temp+(1|source),gaussian,data=data[data$F1parasite=="infection",])
Anova(model,type=3)

model=glmer(ageatfirst~F0temp+F1temp+F0parasite+F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
Anova(model,type=3)

model=glmer(firstclutch~F0temp+F1temp+F0parasite+F1parasite+(1|source),poisson,data=data2[data2$infection.status!="exposed-uninfected",])
Anova(model,type=3)

####parasite traits for initial infection process####
model=glmer(log(Hemocytes.per.spore+1)~F0temp*F1temp*F0parasite+(1|source),gaussian,data=data[data$F1parasite=="infection",])
Anova(model,type=3)

####host traits for terminal infection process####
model=glmer.nb(totalbaby~F0temp*F1temp*F0parasite*F1parasite+(1|source),data=data2[data2$infection.status!="exposed-uninfected",])
Anova(model,type=3)

#considering immune responses in predicting fecundity
model=glmer.nb(totalbaby~F0parasite*(F1temp+F0temp)*(Hemocytes.per.spore)+F1temp*F0temp+(Hemocytes.per.spore)*(F1temp+F0temp)+(gut.resistance)+F1temp*(F0temp+F0parasite)+gut.resistance+(1|source),data=data2[data2$infection=="1",])
Anova(model,type=3)

model=coxme(Surv(lifespan)~F0temp*F1temp*F0parasite*F1parasite+(1|source),data=data2[data2$infection.status!="exposed-uninfected",])
Anova(model,type=3)

####parasite traits for terminal infection process####
model=glmer(infection~F0temp+F1temp+F0parasite+(1|source),binomial,data=data2)
Anova(model,type=3)

model=glmer(log(maturetrans+1)~F0temp+F1temp+F0parasite+gut.resistance+Hemocytes.per.spore+(1|source),gaussian,data=data2[data2$infection=="1",])
Anova(model,type=3)

