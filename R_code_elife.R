############################################################
#T.Constant 
#workspace prep
############################################################

#set directory
setwd(" ")

###load packages

library(caper)
library(apTreeshape)
library(phytools)
library(ape)
library(caper)
library(sjPlot)
library(devtools)
library(nlme)
library(geiger)
library(AICcmodavg)
library(sciplot)
library(MuMIn)
library(car)
library(ggplot2)

#Phylogenetic data

Tree_m_1_2<-read.nexus("model_1_and_2.nex")
Tree_m_3<-read.nexus("model_3.nex")
Tree_m_4<-read.nexus("model_4.nex")
Tree_m_5<-read.nexus("model_5.nex")
Tree_m_6<-read.nexus("model_6.nex")
Tree_m_7<-read.nexus("model_7.nex")


#Biological data 
#all quantitative variables were standardized (using z-scores) in multi-factor models

model_1_and_2<-read.csv2("model_1_and_2.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)
model_3<-read.csv2("model_3.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)

model_4<-read.csv2("model_4.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)
model_5<-read.csv2("model_5.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)
#standardization 
model_5$zlate_mating=scale(model_5$late_mating,center=TRUE,scale = TRUE)
model_5$zlog_rel_testes_mass=scale(model_5$log_rel_testes_mass,center=TRUE,scale = TRUE)
model_5$zmin_temper=scale(model_5$min_temper,center=TRUE,scale=TRUE)

model_6<-read.csv2("model_6.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)
#standardization 
model_6$zlate_mating=scale(model_6$late_mating,center=TRUE,scale = TRUE)
model_6$zbody_mass_before_mating=scale(model_6$body_mass_before_mating,center=TRUE,scale = TRUE)
model_6$zmin_temper=scale(model_6$min_temper,center=TRUE,scale=TRUE)

model_7<-read.csv2("model_7.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)
#standardization 
model_7$ztime_after_mating=scale(model_7$time_after_mating,center=TRUE,scale = TRUE)
model_7$zmaternal_effort=scale(model_7$maternal_effort,center=TRUE,scale = TRUE)


#Build consensus tree from multiphylo object
Tree_m_1_2_cons<-consensus(Tree_m_1_2,p=1,check.labels = TRUE, rooted = TRUE)
Tree_m_1_2_cons

Tree_m_3_cons<-consensus(Tree_m_3,p=1,check.labels = TRUE, rooted = TRUE)
Tree_m_3_cons

Tree_m_4_cons<-consensus(Tree_m_4,p=1,check.labels = TRUE, rooted = TRUE)
Tree_m_4_cons

Tree_m_5_cons<-consensus(Tree_m_5,p=1,check.labels = TRUE, rooted = TRUE)
Tree_m_5_cons

Tree_m_6_cons<-consensus(Tree_m_6,p=1,check.labels = TRUE, rooted = TRUE)
Tree_m_6_cons

Tree_m_7_cons<-consensus(Tree_m_7,p=1,check.labels = TRUE, rooted = TRUE)
Tree_m_7_cons

#Figure for paper

par(mfrow = c(1, 1))
plot(Tree_m_1_2_cons,cex=1,no.margin = TRUE,label.offset = 1)
plot(Tree_m_3_cons,cex=1,no.margin = TRUE,label.offset = 1)
plot(Tree_m_4_cons,cex=1,no.margin = TRUE,label.offset = 1)
plot(Tree_m_5_cons,cex=1,no.margin = TRUE,label.offset = 1)
plot(Tree_m_6_cons,cex=1,no.margin = TRUE,label.offset = 1)
plot(Tree_m_7_cons,cex=1,no.margin = TRUE,label.offset = 1)


#Determine branch lengths of consensus trees 
#Are the tree ultrametric?  (the answer must be TRUE)

Tree_m_1_2_comp <- compute.brlen(Tree_m_1_2_cons, method="Grafen") 
is.ultrametric(Tree_m_1_2_comp)

Tree_m_3_comp <- compute.brlen(Tree_m_3_cons, method="Grafen") 
is.ultrametric(Tree_m_3_comp)

Tree_m_4_comp <- compute.brlen(Tree_m_4_cons, method="Grafen") 
is.ultrametric(Tree_m_4_comp)

Tree_m_5_comp <- compute.brlen(Tree_m_5_cons, method="Grafen") 
is.ultrametric(Tree_m_5_comp)

Tree_m_6_comp <- compute.brlen(Tree_m_6_cons, method="Grafen") 
is.ultrametric(Tree_m_6_comp)

Tree_m_7_comp <- compute.brlen(Tree_m_7_cons, method="Grafen") 
is.ultrametric(Tree_m_7_comp)


#Are the tree dichotomus? (the answer must be TRUE)

is.binary.multiPhylo(Tree_m_1_2)
is.binary.tree(Tree_m_1_2_comp)

is.binary.multiPhylo(Tree_m_3)
is.binary.tree(Tree_m_3_comp)

is.binary.multiPhylo(Tree_m_4)
is.binary.tree(Tree_m_4_comp)

is.binary.multiPhylo(Tree_m_5)
is.binary.tree(Tree_m_5_comp)

is.binary.multiPhylo(Tree_m_6)
is.binary.tree(Tree_m_6_comp)

is.binary.multiPhylo(Tree_m_7)
is.binary.tree(Tree_m_7_comp)

#Combine both datafiles (phylo and Biological data)

Tree_m_1_2_comp$node.label<-NULL
comb_model_1_and_2<-comparative.data(phy=Tree_m_1_2_comp,data=model_1_and_2,names.col=species,vcv=TRUE,na.omit=FALSE,warn.dropped=TRUE)
str(comb_model_1_and_2)
sort(unique(Tree_m_1_2_comp$tip.label))
sort(unique(model_1_and_2$species))
sort(unique(comb_model_1_and_2$phy$tip.label))
comb_model_1_and_2$dropped$tips
sort(unique(comb_model_1_and_2$phy$tip.label))==sort(unique(model_1_and_2$species))

Tree_m_3_comp$node.label<-NULL
comb_model_3<-comparative.data(phy=Tree_m_3_comp,data=model_3,names.col=species,vcv=TRUE,na.omit=FALSE,warn.dropped=TRUE)
str(comb_model_3)
sort(unique(Tree_m_3_comp$tip.label))
sort(unique(model_3$species))
sort(unique(comb_model_3$phy$tip.label))
comb_model_3$dropped$tips
sort(unique(comb_model_3$phy$tip.label))==sort(unique(model_3$species))

Tree_m_4_comp$node.label<-NULL
comb_model_4<-comparative.data(phy=Tree_m_4_comp,data=model_4,names.col=species,vcv=TRUE,na.omit=FALSE,warn.dropped=TRUE)
str(comb_model_4)
sort(unique(Tree_m_4_comp$tip.label))
sort(unique(model_4$species))
sort(unique(comb_model_4$phy$tip.label))
comb_model_4$dropped$tips
sort(unique(comb_model_4$phy$tip.label))==sort(unique(model_4$species))

Tree_m_5_comp$node.label<-NULL
comb_model_5<-comparative.data(phy=Tree_m_5_comp,data=model_5,names.col=species,vcv=TRUE,na.omit=FALSE,warn.dropped=TRUE)
str(comb_model_5)
sort(unique(Tree_m_5_comp$tip.label))
sort(unique(model_5$species))
sort(unique(comb_model_5$phy$tip.label))
comb_model_5$dropped$tips
sort(unique(comb_model_5$phy$tip.label))==sort(unique(model_5$species))

Tree_m_6_comp$node.label<-NULL
comb_model_6<-comparative.data(phy=Tree_m_6_comp,data=model_6,names.col=species,vcv=TRUE,na.omit=FALSE,warn.dropped=TRUE)
str(comb_model_6)
sort(unique(Tree_m_6_comp$tip.label))
sort(unique(model_6$species))
sort(unique(comb_model_6$phy$tip.label))
comb_model_6$dropped$tips
sort(unique(comb_model_6$phy$tip.label))==sort(unique(model_6$species))

Tree_m_7_comp$node.label<-NULL
comb_model_7<-comparative.data(phy=Tree_m_7_comp,data=model_7,names.col=species,vcv=TRUE,na.omit=FALSE,warn.dropped=TRUE)
str(comb_model_7)
sort(unique(Tree_m_7_comp$tip.label))
sort(unique(model_7$species))
sort(unique(comb_model_7$phy$tip.label))
comb_model_7$dropped$tips
sort(unique(comb_model_7$phy$tip.label))==sort(unique(model_7$species))

#PGLS models

#pgls model with lambda evaluation by ML (1=high phylogenetic covariance; 0 = no phylogenetic covariance)

#model 1

fit1<-pgls(body_mass_before_mating~body_mass_during_mating
           
           ,data=comb_model_1_and_2,lambda = "ML")

summary(fit1)

#Normality and homoscedasticity are checked by graphical observation

#1) normality
#Ploting normal qqplot
#If the data are normally distributed, the points on a Q-Q graph 
#will lie on a straight diagonal line.

#2) homoscedasticity 
#Ploting residual value vs fitted value allows to check the homoscedasticity
#The residuals have to be randomly distributed around the 0 line, 
#approximately form a horizontal band around the 0 line 
#with no extreme value from the basic random pattern of residuals. 

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit1)


#model 2

fit2<-pgls(body_mass_before_mating~log_rel_testes_mass
                      
          ,data=comb_model_1_and_2,lambda = "ML")
           
summary(fit2)

#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit2)

#model 3

fit3<-pgls(time_after_mating~body_mass_during_mating
           
           ,data=comb_model_3,lambda = "ML")

summary(fit3)

#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit3)


#model 3 without phylogenetic covariance

fit3a<-pgls(time_after_mating~body_mass_during_mating
            
            ,data=comb_model_3,lambda = 0.00001)

summary(fit3a)

#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit3a)


#model 4

fit4<-pgls(time_after_mating~body_mass_before_during_mating
           
           ,data=comb_model_4,lambda = "ML")

summary(fit4)

#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit4)

#model 4 without phylogenetic covariance

fit4a<-pgls(time_after_mating~body_mass_before_during_mating
            
            ,data=comb_model_4,lambda = 0.00001)

summary(fit4a)

#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit4a)

#model 5
#We test the full model and keep the best model based on the AICc

fit5<-pgls(protandry~
             zlog_rel_testes_mass*min_temper+
             zlog_rel_testes_mass*foodstoring+
             zlog_rel_testes_mass*late_mating
           ,data=comb_model_5,lambda = "ML")

summary(model.avg(dredge(fit5, rank="AICc"),delta<6))

#best model 5 :

fit5<-pgls(protandry~
             zlog_rel_testes_mass*zmin_temper
           ,data=comb_model_5,lambda = "ML")

summary(fit5)

#We tested multicollinearity using the variance inflation factor (VIF < 3) 
#on linear models including the factors of the best models

fit5vif<-lm(protandry~
             zlog_rel_testes_mass*zmin_temper
           ,data=model_5)
vif(fit5vif)

#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit5)


#graph interaction with minimum temperature and without standardization for a better illustration

fit5<-pgls(protandry~
             log_rel_testes_mass*min_temper
           ,data=comb_model_5,lambda = "ML")

summary(fit5)

q1<-min(model_5$min_temper)
q1
q2<-mean(model_5$min_temper)
q2
q3<-max(model_5$min_temper)
q3

eq4 = function(x){ 25.22926  +(1.01583   *-29.85)+(((3.02791)+(1.50432*-29.85))*x)}
eq2 = function(x){25.22926  +(1.01583  *-8.971818)+(((3.02791   )+(1.50432*-8.971818))*x)}
eq = function(x){25.22926 +(1.01583  *14.6)+(((3.02791   )+(1.50432*14.6))*x)}
mid<-mean(model_5$min_temper)
ggplot(model_5, aes(x=log_rel_testes_mass, y=protandry,colour=min_temper))+
  geom_point(size=4)+
  theme_classic()+
  scale_color_gradient2("Min temperature (°C)",midpoint=mid,  low="blue4", mid="snow3",
                        high="red4",  space = "Lab")+
  geom_function(fun=eq, color="red4")+
  geom_function(fun=eq4, color="blue4")+
  geom_function(fun=eq2, color="snow3")+
  ylab("Protandry (day)")+
  xlab("Log(Relative testes mass)")+
  ylim(0,50)


#best model 5 without phylogenetic covariance


fit5a<-pgls(protandry~zlog_rel_testes_mass*zmin_temper
          ,data=comb_model_5,lambda = 0.00001)

summary(fit5a)

#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit5a)

#model 6
#We test the full model and keep the best model based on the AICc

fit6<-pgls(protandry~
             zbody_mass_before_mating*min_temper+
             zbody_mass_before_mating*foodstoring+
             zbody_mass_before_mating*late_mating
           ,data=comb_model_6,lambda = "ML")

summary(model.avg(dredge(fit6, rank="AICc"),delta<6))

#best model 6 :

fit6<-pgls(protandry~
             zbody_mass_before_mating+zlate_mating
           ,data=comb_model_6,lambda = "ML")

summary(fit6)

#We tested multicollinearity using the variance inflation factor (VIF < 3) 
#on linear models including the factors of the best models

fit6vif<-lm(protandry~
              zbody_mass_before_mating+zlate_mating
            ,data=model_6)
vif(fit6vif)


#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit6)

#graph model 6

ggplot(model_6, aes(x=zbody_mass_before_mating, y=protandry,colour=late_mating))+
  geom_point(size=4)+
  theme_classic()+
  
  
  labs(x="Standardized-Body mass change before mating",
       y="Protandry (day)")+
  scale_y_continuous(limits = c(-60,10), expand = c(0, 0)) +
  scale_color_gradient("Late mating (week)")+
  geom_abline(intercept = -24.9654 , slope =-5.9514,color="red", 
              size=1)

#graph model 6 without standardization for a better illustration

fit6b<-pgls(protandry~
             body_mass_before_mating+late_mating
           ,data=comb_model_6,lambda = "ML")

summary(fit6b)

#determine the average value of late mating
q4<-mean(model_6$late_mating)
q4

#equation for regression line of model 6 with mean value of late mating

eq4 = function(x){ 31.77956  +(-4.89331   *1.815789)+(0.54228*x)}

ggplot(model_6, aes(x=body_mass_before_mating, y=protandry,colour=late_mating))+
  geom_point(size=4)+
  theme_classic()+
  
  
  labs(x="Body mass change before mating (%)",
       y="Protandry (day)")+
  scale_y_continuous(limits = c(0,55), expand = c(0, 0)) +
  scale_color_gradient("Late mating (week)")+
  geom_function(fun=eq4, color="red",size=1)


#best model 6 without phylogenetic covariance


fit6a<-pgls(protandry~ zbody_mass_before_mating+zlate_mating
            ,data=comb_model_6,lambda = 0.00001)

summary(fit6a)

#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit6a)

#model 7
#We test the full model and keep the best model based on the AICc

fit7<-pgls(sex_diff_immerg~
             time_after_mating+maternal_effort
           ,data=comb_model_7,lambda = "ML")

summary(model.avg(dredge(fit7, rank="AICc")))

#best model 7 :

fit7<-pgls(sex_diff_immerg~
             ztime_after_mating+zmaternal_effort
           ,data=comb_model_7,lambda = "ML")

summary(fit7)


#We tested multicollinearity using the variance inflation factor (VIF < 3) 
#on linear models including the factors of the best models

fit7vif<-lm(sex_diff_immerg~
              ztime_after_mating+zmaternal_effort
            ,data=model_7)
vif(fit7vif)


#graph model 7 without standardization for a better illustration

fit7b<-pgls(sex_diff_immerg~
             time_after_mating+maternal_effort
           ,data=comb_model_7,lambda = "ML")

summary(fit7b)

#determine the average value of maternal effort
q5<-mean(model_7$maternal_effort)
q5

#equation for regression line of model 6 with mean value of late mating

eq5 = function(x){ 12.433706  +(0.420307   *67.1825)+(-0.344105*x)}

ggplot(model_7, aes(x=time_after_mating, y=sex_diff_immerg,colour=maternal_effort))+
  geom_point(size=4)+
  theme_classic()+
  scale_color_gradientn("Maternal effort (day)",colours = rainbow(3))+
  ylab("Sex difference in immergence (day)")+
  xlab("Post-mating activity time (day)")+
  geom_function(fun=eq5, color="red",size=1)
