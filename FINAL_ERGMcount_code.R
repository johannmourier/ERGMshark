rm(list=ls())

################
#######Requiered libraries
library(sna)
library(statnet)
library (intergraph)
library(igraph)
library(gdistance)
library (fastnet)
library(ergm)
library(ergm.count)
library(stringr)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(colorRamps)
library(prettymapr)
library(wesanderson)
library(plotfunctions)
library(party)
library(randomForest)
library(randomForestExplainer)

#Added a couple of new packages
library(coda)
library(asbio)

# Set the working directory 
setwd("/Users/johannmourier/Documents/Scripts/CHARC/Stage_Angelique/Etape 7 - Modele ERGM/NEW")


######################################################################################
############# STEP 1 ##################################################
#######################################################################

############# Prepare parameters data as input for models

param=get(load("Decision param-V3.RData")) #Array with: stations / parameters / months

requins=get(load("Decision requins.RData")) # Data frame with: sharks / months
stations=get(load("Decision stations-V3.RData")) # Data frame with: stations / months

mouvments=get(load("mouvments-V3.RData")) # Array withshark movements : stations / stations / sharks / months

Geodist=get(load("Geodist-withDCP.RData")) # Matrix of geographic distance between stations : stations / stations
diag(Geodist)<-0

# Load dataset
bdd=get(load("Baseded-V3-sharkpres.RData")) 

tableIbis=get(load("Baseded-V3-Envir.RData"))

requins<-requins[,c(12:41)] #Remove useless months

stations<-stations[dimnames(param)[1][[1]],] #Order stations following a same order

# Get names of sharks and stations
shark.names <- rownames(requins)
station.names <- rownames(stations)
mois.names <- colnames(requins)

# 39 sharks in shark.names and only 32 in the database, because residency is only available 
# for 32 sharks. We restrict the analysis for 32 sharks
names<-read.csv("shark.name.csv",sep=",",header=T)
bdd<-merge(bdd,names, by=c("nom_requin","code_esp"))

shark.names<- shark.names[which(shark.names %in% as.character(unique(bdd$nom_requin)))]


#### create matrix : months / stations / sharks / residency with same-different sexe
diffsexe<-array(0,dim=c(length(mois.names),length(station.names),length(shark.names),2),dimnames=list(mois.names,station.names,shark.names,c("TP_Meme","TP_Diff")))
for (i in shark.names) { #Choice of focal shark
  for (j in mois.names){
    sexe.focal<- as.character(unique(bdd[ bdd$nom_requin==i,"sex"])) #sex of focal shark
    tp.m<-param[,40,j]
    tp.f<-param[,39,j]
    if( sexe.focal == "M"){TP_Meme<- tp.m ;  TP_Diff<- tp.f} #Definiton of residency (TP=presence time) of same (TP Meme) or different (TP Diff) sex according to the sex of the forcal shark  
    if( sexe.focal == "F"){ TP_Diff<- tp.m ;  TP_Meme<- tp.f}
    if(length(TP_Meme)>0) {diffsexe[j,,i,1]<-TP_Meme} #Storage in the array diffsexe
    if(length(TP_Diff)>0) {diffsexe[j,,i,2]<-TP_Diff}
  }
}


#Names of stations, parameters, months 
station.nom<-dimnames(param)[[1]]
param.nom<-dimnames(param)[[2]]
mois.nom<-unlist(dimnames(param)[3])

dim1<-length(station.nom)
dim2<-length(param.nom) 
dim3<-length(mois.nom)
param1<-array(NA,dim=c(dim1,dim2,dim3),dimnames=list(station.nom,param.nom,mois.nom))

# Scale environnemental data
for (p in param.nom){
  for (j in as.numeric(mois.nom)){
    for (s in station.nom){
      param1[s,p,j]<-(param[s,p,j]-mean(param[,p,],na.rm=T))/(sd(param[,p,],na.rm=T))
    }
  }
}

save(param1,file="Param scale3.RData")

dim2<-length(station.nom)
dim3<-length(shark.names) 
dim1<-length(mois.nom)
diff.sex<-array(0,dim=c(length(mois.names),length(station.names),length(shark.names),2),dimnames=list(mois.names,station.names,shark.names,c("TP_Meme","TP_Diff")))
for (i in as.numeric(mois.nom)){
  for (j in station.nom){
    for (s in shark.names){
      diff.sex[i,j,s,1]<-(diffsexe[i,j,s,1]-mean(diffsexe[,,s,1],na.rm=T))/(sd(diffsexe[,,s,1],na.rm=T))
      diff.sex[i,j,s,2]<-(diffsexe[i,j,s,2]-mean(diffsexe[,,s,2],na.rm=T))/(sd(diffsexe[,,s,2],na.rm=T))
    }}}

data_param=get(load("Param scale3.RData"))




######################################################################################
############# STEP 2 ##################################################
#######################################################################

### Visialize environemental data variability across locations
library(reshape)
# reorganize the table

#Rain (meanRR.mm.)
rain<-as.data.frame(param[,"meanRR.mm.",])
rain<-cbind(rownames(param),rain)
Rain<-melt(rain)
colnames(Rain)<-c("Receiver","Months","Rain")
Rain$Data<-"Rain"

g1<-ggplot(Rain, aes(x=Months, y=Rain)) +
  geom_boxplot(colour = "gray50",fill = "gray70", alpha = 0.3)+
  geom_point(colour = "#FC4E07", size=1)+
  theme_bw() 

# Swell (MeanHsig)
swell<-as.data.frame(param[,"MeanHsig",])
swell<-cbind(rownames(param),swell)
Swell<-melt(swell)
colnames(Swell)<-c("Receiver","Months","Swell")
Swell$Data<-"Swell"

g2<-ggplot(Swell, aes(x=Months, y=Swell)) +
  geom_boxplot(colour = "gray50",fill = "gray70", alpha = 0.3)+
  geom_point(colour = "#FC4E07", size=1)+
  theme_bw() 

# Mean reflectance (rs_moy)
rs<-as.data.frame(param[,"rs_moy",])
rs<-cbind(rownames(param),rs)
Reflectance<-melt(rs)
colnames(Reflectance)<-c("Receiver","Months","Reflectance")
Reflectance$Data<-"Reflectance"

g3<-ggplot(Reflectance, aes(x=Months, y=Reflectance)) +
  geom_boxplot(colour = "gray50",fill = "gray70", alpha = 0.3)+
  geom_point(colour = "#FC4E07", size=1)+
  theme_bw() 

# Turtle density (MeanDensityTurtle)
turtle<-as.data.frame(param[,"MeanDensityTurtle",])
turtle<-cbind(rownames(param),turtle)
Turtle<-melt(turtle)
colnames(Turtle)<-c("Receiver","Months","Turtle")
Turtle$Data<-"Turtle"

g4<-ggplot(Turtle, aes(x=Months, y=Turtle)) +
  geom_boxplot(colour = "gray50",fill = "gray70", alpha = 0.3)+
  geom_point(colour = "#FC4E07", size=1)+
  theme_bw() 

MeanTotalFreq

# Human Frequentation (MeanTotalFreq)
anthro<-as.data.frame(param[,"MeanTotalFreq",])
anthro<-cbind(rownames(param),anthro)
Anthro<-melt(anthro)
colnames(Anthro)<-c("Receiver","Months","Human_use")
Anthro$Data<-"Human_use"

g5<-ggplot(Anthro, aes(x=Months, y=Human_use)) +
  geom_boxplot(colour = "gray50",fill = "gray70", alpha = 0.3)+
  geom_point(colour = "#FC4E07", size=1)+
  theme_bw() 

# Tiger shark occupancy (TxoccTS)
TSocc<-as.data.frame(param[,"TxoccTS",])
TSocc<-cbind(rownames(param),TSocc)
TSocc<-melt(TSocc)
colnames(TSocc)<-c("Receiver","Months","TigerOccupancy")
TSocc$Data<-"TigerOccupancy"

g6<-ggplot(TSocc, aes(x=Months, y=TigerOccupancy)) +
  geom_boxplot(colour = "gray50",fill = "gray70", alpha = 0.3)+
  geom_point(colour = "#FC4E07", size=1)+
  theme_bw() 

# Bull shark occupancy (TxoccBS)
BSocc<-as.data.frame(param[,"TxoccBS",])
BSocc<-cbind(rownames(param),BSocc)
BSocc<-melt(BSocc)
colnames(BSocc)<-c("Receiver","Months","BullOccupancy")
BSocc$Data<-"BullOccupancy"

g7<-ggplot(BSocc, aes(x=Months, y=BullOccupancy)) +
  geom_boxplot(colour = "gray50",fill = "gray70", alpha = 0.3)+
  geom_point(colour = "#FC4E07", size=1)+
  theme_bw() 

library(ggpubr)
pdf(file = "Variables.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 12) # The height of the plot in inches

figure <- ggarrange(g1, g2, g3, g4, g5, g6, g7,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 7)
figure
dev.off()



######################################################################################
############# STEP 3 ##################################################
#######################################################################

#########################     NETWORK SELECTION      ##########################

time.cutoff<-20 # The shark needs to have an active tag during at least 20 days within the month (e.g can't be tagged at the end of the month)
station.cutoff<-20 # Receivers (stations) need to be deployed at least on 20 days during a month 

# Create a table with ERGM count outputs 
var<-c("TxoccBS","TxoccTS","MeanHsig","MeanDensityTurtle","meanRR.mm.","MeanTotalFreq",
       "rs_moy","nb_bull","nb_tiger","TP_Meme","TP_Diff")

l<- 1
##### Create a list of network objects
networks<-list()
names<-list()
months<-list()
####### Procedure that select the terms of the ergm model for each shark and each month
for (i in shark.names ){
  date<-colnames(requins)[which(requins[i,]>time.cutoff)] # month in which the shark is tagged at least 20 days
  if (length(date)>0){
    for (j in date){
      stat<-rownames(stations)[which(as.numeric(as.character(stations[,j]))>time.cutoff)] #stations deployees plus de 20 jours durant le mois j
      env <-param1[stat,,j] # Parameter measured in these stations during month j 
      env<-env[,c("TxoccBS","TxoccTS","meanRR.mm.","MeanHsig","MeanDensityTurtle","MeanTotalFreq",
                  "rs_moy","nb_bull","nb_tiger")] 
      if(length(stat)>0){
        mov.tmp<-mouvments[stat,stat,i,j]
        env<-na.omit(env)
        mov.tmp<-mov.tmp[ rownames(env), rownames(env)]  
        if(nrow(env)>station.cutoff){ # verify that more than 20 stations are taken in the network 
          
          #Create a network containing movements of shark i in month j 
          g<-as.network(mov.tmp, matrix.type = "adjacency", directed = TRUE,names.eval = "weight",ignore.eval = FALSE,loops=F)
          ver.attr<-(rownames(env))[which((rownames(env)) %in% (g %v% "vertex.names"))] # Harmonization of station selection
          ver.attr2<-(rownames(env))[which((rownames(env)) %in% (dimnames(diff.sex)[2][[1]]))]
          ver.attr1<-env[ver.attr,] #Contains parameter measurements for selected stations
          ver.attr1<-cbind(ver.attr1,TP_Meme=diff.sex[j,ver.attr2,i,1]) # Add residency CRT of same and different sex of focal shark 
          ver.attr1<-cbind(ver.attr1,TP_Diff=diff.sex[j,ver.attr2,i,2])
          ver.attr1[is.na(ver.attr1)]<-0
          # Add parameters as node attributes 
          g %v% "TxoccBS" <- ver.attr1[,1]  
          g %v% "TxoccTS" <- ver.attr1[,2]
          g %v% "meanRR.mm." <- ver.attr1[,3]
          g %v% "MeanHsig" <- ver.attr1[,4]
          g %v% "MeanDensityTurtle" <- ver.attr1[,5]
          g %v% "MeanTotalFreq"<- ver.attr1[,6]
          g %v% "rs_moy"<- ver.attr1[,7]
          g %v% "nb_bull"<- ver.attr1[,8]
          g %v% "nb_tiger"<- ver.attr1[,9]
          g %v% "TP_Meme"<-ver.attr1[,10]
          g %v% "TP_Diff"<-ver.attr1[,11]
          g %e% "Geodist"<-Geodist[rownames(ver.attr1),rownames(ver.attr1)] 
          ## we keep more than 5 edges
          if(network.edgecount(g) >=10){
            
            networks[[l]]<-g 
            names[[l]]<-i      
            months[[l]]<-j      
            mylist <- list( networks , names , months)
            l<- l+1
            rm(list = c("stat","env","g","ver.attr","ver.attr1"))
            
          }
        }
        
      } 
    }
    print(paste("Le requin",i,"mois",j))
  }
  rm(list="date")
}



######################################################################################
############# STEP 4 ##################################################
#######################################################################

#################### MODEL SELECTION ########################
#############################################################

### for each network we select the best model from smallest AIC value

####
### for term 'nodeicov' added to the list of the networks (6th dimension)
for (i in 1:length(networks)){
  g<-networks[[i]]
  t_AICs<-rep(NA,68)
  #options(warn=2)
  x1<-try(ergm(g~nonzero+sum+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x1)=="ergm")){t_AICs[1]<-AIC(x1)}
  x2<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x2)!="try-error")){t_AICs[2]<-AIC(x2)}
  x3<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x3)!="try-error")){t_AICs[3]<-AIC(x3)}
  x4<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+nodeicov('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x4)!="try-error")){t_AICs[4]<-AIC(x4)}
  x5<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x5)!="try-error")){t_AICs[5]<-AIC(x5)}
  x6<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x6)!="try-error")){t_AICs[6]<-AIC(x6)}
  x7<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x7)!="try-error")){t_AICs[7]<-AIC(x7)}
  x8<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x8)!="try-error")){t_AICs[8]<-AIC(x8)}
  x9<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x9)!="try-error")){t_AICs[9]<-AIC(x9)}
  x10<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x10)!="try-error")){t_AICs[10]<-AIC(x10)}
  x11<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x11)!="try-error")){t_AICs[11]<-AIC(x11)}
  x12<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x12)!="try-error")){t_AICs[12]<-AIC(x12)}
  x13<-try(ergm(g~nonzero+sum+nodeicov('TxoccBS')+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x13)!="try-error")){t_AICs[13]<-AIC(x13)}
  x14<-try(ergm(g~nonzero+sum+nodeicov('TxoccTS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x14)!="try-error")){t_AICs[14]<-AIC(x14)}
  x15<-try(ergm(g~nonzero+sum+nodeicov('TxoccTS')+nodeicov('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x15)!="try-error")){t_AICs[15]<-AIC(x15)}
  x16<-try(ergm(g~nonzero+sum+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x16)!="try-error")){t_AICs[16]<-AIC(x16)}
  x17<-try(ergm(g~nonzero+sum+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x17)!="try-error")){t_AICs[17]<-AIC(x17)}
  x18<-try(ergm(g~nonzero+sum+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x18)!="try-error")){t_AICs[18]<-AIC(x18)}
  x19<-try(ergm(g~nonzero+sum+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x19)!="try-error")){t_AICs[19]<-AIC(x19)}
  x20<-try(ergm(g~nonzero+sum+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x20)!="try-error")){t_AICs[20]<-AIC(x20)}
  x21<-try(ergm(g~nonzero+sum+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x21)!="try-error")){t_AICs[21]<-AIC(x21)}
  x22<-try(ergm(g~nonzero+sum+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x22)!="try-error")){t_AICs[22]<-AIC(x22)}
  x23<-try(ergm(g~nonzero+sum+nodeicov('TxoccTS')+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x23)!="try-error")){t_AICs[23]<-AIC(x23)}
  x24<-try(ergm(g~nonzero+sum+nodeicov('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x24)!="try-error")){t_AICs[24]<-AIC(x24)}
  x25<-try(ergm(g~nonzero+sum+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x25)!="try-error")){t_AICs[25]<-AIC(x25)}
  x26<-try(ergm(g~nonzero+sum+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x26)!="try-error")){t_AICs[26]<-AIC(x26)}
  x27<-try(ergm(g~nonzero+sum+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x27)!="try-error")){t_AICs[27]<-AIC(x27)}
  x28<-try(ergm(g~nonzero+sum+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x28)!="try-error")){t_AICs[28]<-AIC(x28)}
  x29<-try(ergm(g~nonzero+sum+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x29)!="try-error")){t_AICs[29]<-AIC(x29)}
  x30<-try(ergm(g~nonzero+sum+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x30)!="try-error")){t_AICs[30]<-AIC(x30)}
  x31<-try(ergm(g~nonzero+sum+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x31)!="try-error")){t_AICs[31]<-AIC(x31)}
  x32<-try(ergm(g~nonzero+sum+nodeicov('MeanHsig')+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x32)!="try-error")){t_AICs[32]<-AIC(x32)}
  x33<-try(ergm(g~nonzero+sum+nodeicov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x33)!="try-error")){t_AICs[33]<-AIC(x33)}
  x34<-try(ergm(g~nonzero+sum+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x34)!="try-error")){t_AICs[34]<-AIC(x34)}
  x35<-try(ergm(g~nonzero+sum+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x35)!="try-error")){t_AICs[35]<-AIC(x35)}
  x36<-try(ergm(g~nonzero+sum+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x36)!="try-error")){t_AICs[36]<-AIC(x36)}
  x37<-try(ergm(g~nonzero+sum+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x37)!="try-error")){t_AICs[37]<-AIC(x37)}
  x38<-try(ergm(g~nonzero+sum+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x38)!="try-error")){t_AICs[38]<-AIC(x38)}
  x39<-try(ergm(g~nonzero+sum+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x39)!="try-error")){t_AICs[39]<-AIC(x39)}
  x40<-try(ergm(g~nonzero+sum+nodeicov('meanRR.mm.')+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x40)!="try-error")){t_AICs[40]<-AIC(x40)}
  x41<-try(ergm(g~nonzero+sum+nodeicov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x41)!="try-error")){t_AICs[41]<-AIC(x41)}
  x42<-try(ergm(g~nonzero+sum+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x42)!="try-error")){t_AICs[42]<-AIC(x42)}
  x43<-try(ergm(g~nonzero+sum+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x43)!="try-error")){t_AICs[43]<-AIC(x43)}
  x44<-try(ergm(g~nonzero+sum+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x44)!="try-error")){t_AICs[44]<-AIC(x44)}
  x45<-try(ergm(g~nonzero+sum+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x45)!="try-error")){t_AICs[45]<-AIC(x45)}
  x46<-try(ergm(g~nonzero+sum+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x46)!="try-error")){t_AICs[46]<-AIC(x46)}
  x47<-try(ergm(g~nonzero+sum+nodeicov('rs_moy')+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x47)!="try-error")){t_AICs[47]<-AIC(x47)}
  x48<-try(ergm(g~nonzero+sum+nodeicov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x48)!="try-error")){t_AICs[48]<-AIC(x48)}
  x49<-try(ergm(g~nonzero+sum+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x49)!="try-error")){t_AICs[49]<-AIC(x49)}
  x50<-try(ergm(g~nonzero+sum+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x50)!="try-error")){t_AICs[50]<-AIC(x50)}
  x51<-try(ergm(g~nonzero+sum+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x51)!="try-error")){t_AICs[51]<-AIC(x51)}
  x52<-try(ergm(g~nonzero+sum+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x52)!="try-error")){t_AICs[52]<-AIC(x52)}
  x53<-try(ergm(g~nonzero+sum+nodeicov('MeanDensityTurtle')+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x53)!="try-error")){t_AICs[53]<-AIC(x53)}
  x54<-try(ergm(g~nonzero+sum+nodeicov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x54)!="try-error")){t_AICs[54]<-AIC(x54)}
  x55<-try(ergm(g~nonzero+sum+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x55)!="try-error")){t_AICs[55]<-AIC(x55)}
  x56<-try(ergm(g~nonzero+sum+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x56)!="try-error")){t_AICs[56]<-AIC(x56)}
  x57<-try(ergm(g~nonzero+sum+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x57)!="try-error")){t_AICs[57]<-AIC(x57)}
  x58<-try(ergm(g~nonzero+sum+nodeicov('MeanTotalFreq')+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x58)!="try-error")){t_AICs[58]<-AIC(x58)}
  x59<-try(ergm(g~nonzero+sum+nodeicov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x59)!="try-error")){t_AICs[59]<-AIC(x59)}
  x60<-try(ergm(g~nonzero+sum+nodeicov('nb_bull')+nodeicov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x60)!="try-error")){t_AICs[60]<-AIC(x60)}
  x61<-try(ergm(g~nonzero+sum+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x61)!="try-error")){t_AICs[61]<-AIC(x61)}
  x62<-try(ergm(g~nonzero+sum+nodeicov('nb_bull')+nodeicov('nb_tiger')+nodeicov('TP_Meme')+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x62)!="try-error")){t_AICs[62]<-AIC(x62)}
  x63<-try(ergm(g~nonzero+sum+nodeicov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x63)!="try-error")){t_AICs[63]<-AIC(x63)}
  x64<-try(ergm(g~nonzero+sum+nodeicov('nb_tiger')+nodeicov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x64)!="try-error")){t_AICs[64]<-AIC(x64)}
  x65<-try(ergm(g~nonzero+sum+nodeicov('nb_tiger')+nodeicov('TP_Meme')+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x65)!="try-error")){t_AICs[65]<-AIC(x65)}
  x66<-try(ergm(g~nonzero+sum+nodeicov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x66)!="try-error")){t_AICs[66]<-AIC(x66)}
  x67<-try(ergm(g~nonzero+sum+nodeicov('TP_Meme')+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x67)!="try-error")){t_AICs[67]<-AIC(x67)}
  x68<-try(ergm(g~nonzero+sum+nodeicov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x68)!="try-error")){t_AICs[68]<-AIC(x68)}
  
  options(warn=1)
  
  b<-which.min(t_AICs)
  
  #b<-which.min(sapply(list(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, #x16,
  #                             x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30#, x31,
  #                             x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45#, x46, x47,
  #                             x48, x49, x50, x51, x52, x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63,
  #                             x64, x65, x66), AIC))
  
  networks[[i]][[6]]<-formula(paste0("x",b))
  
}

# networks[[2]][[6]]

save(networks, file="ergm_nodeicov.RData")

load("ergm_nodeicov.RData")



####
### for term 'nodeocov' added to the list of the networks (7th dimension)

for (i in 1:length(networks)){
  g<-networks[[i]]
  t_AICs<-rep(NA,68)
  #options(warn=2)
  x1<-try(ergm(g~nonzero+sum+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x1)=="ergm")){t_AICs[1]<-AIC(x1)}
  x2<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x2)!="try-error")){t_AICs[2]<-AIC(x2)}
  x3<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x3)!="try-error")){t_AICs[3]<-AIC(x3)}
  x4<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+nodeocov('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x4)!="try-error")){t_AICs[4]<-AIC(x4)}
  x5<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x5)!="try-error")){t_AICs[5]<-AIC(x5)}
  x6<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x6)!="try-error")){t_AICs[6]<-AIC(x6)}
  x7<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x7)!="try-error")){t_AICs[7]<-AIC(x7)}
  x8<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x8)!="try-error")){t_AICs[8]<-AIC(x8)}
  x9<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x9)!="try-error")){t_AICs[9]<-AIC(x9)}
  x10<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x10)!="try-error")){t_AICs[10]<-AIC(x10)}
  x11<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x11)!="try-error")){t_AICs[11]<-AIC(x11)}
  x12<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x12)!="try-error")){t_AICs[12]<-AIC(x12)}
  x13<-try(ergm(g~nonzero+sum+nodeocov('TxoccBS')+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x13)!="try-error")){t_AICs[13]<-AIC(x13)}
  x14<-try(ergm(g~nonzero+sum+nodeocov('TxoccTS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x14)!="try-error")){t_AICs[14]<-AIC(x14)}
  x15<-try(ergm(g~nonzero+sum+nodeocov('TxoccTS')+nodeocov('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x15)!="try-error")){t_AICs[15]<-AIC(x15)}
  x16<-try(ergm(g~nonzero+sum+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x16)!="try-error")){t_AICs[16]<-AIC(x16)}
  x17<-try(ergm(g~nonzero+sum+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x17)!="try-error")){t_AICs[17]<-AIC(x17)}
  x18<-try(ergm(g~nonzero+sum+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x18)!="try-error")){t_AICs[18]<-AIC(x18)}
  x19<-try(ergm(g~nonzero+sum+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x19)!="try-error")){t_AICs[19]<-AIC(x19)}
  x20<-try(ergm(g~nonzero+sum+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x20)!="try-error")){t_AICs[20]<-AIC(x20)}
  x21<-try(ergm(g~nonzero+sum+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x21)!="try-error")){t_AICs[21]<-AIC(x21)}
  x22<-try(ergm(g~nonzero+sum+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x22)!="try-error")){t_AICs[22]<-AIC(x22)}
  x23<-try(ergm(g~nonzero+sum+nodeocov('TxoccTS')+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x23)!="try-error")){t_AICs[23]<-AIC(x23)}
  x24<-try(ergm(g~nonzero+sum+nodeocov('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x24)!="try-error")){t_AICs[24]<-AIC(x24)}
  x25<-try(ergm(g~nonzero+sum+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x25)!="try-error")){t_AICs[25]<-AIC(x25)}
  x26<-try(ergm(g~nonzero+sum+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x26)!="try-error")){t_AICs[26]<-AIC(x26)}
  x27<-try(ergm(g~nonzero+sum+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x27)!="try-error")){t_AICs[27]<-AIC(x27)}
  x28<-try(ergm(g~nonzero+sum+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x28)!="try-error")){t_AICs[28]<-AIC(x28)}
  x29<-try(ergm(g~nonzero+sum+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x29)!="try-error")){t_AICs[29]<-AIC(x29)}
  x30<-try(ergm(g~nonzero+sum+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x30)!="try-error")){t_AICs[30]<-AIC(x30)}
  x31<-try(ergm(g~nonzero+sum+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x31)!="try-error")){t_AICs[31]<-AIC(x31)}
  x32<-try(ergm(g~nonzero+sum+nodeocov('MeanHsig')+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x32)!="try-error")){t_AICs[32]<-AIC(x32)}
  x33<-try(ergm(g~nonzero+sum+nodeocov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x33)!="try-error")){t_AICs[33]<-AIC(x33)}
  x34<-try(ergm(g~nonzero+sum+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x34)!="try-error")){t_AICs[34]<-AIC(x34)}
  x35<-try(ergm(g~nonzero+sum+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x35)!="try-error")){t_AICs[35]<-AIC(x35)}
  x36<-try(ergm(g~nonzero+sum+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x36)!="try-error")){t_AICs[36]<-AIC(x36)}
  x37<-try(ergm(g~nonzero+sum+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x37)!="try-error")){t_AICs[37]<-AIC(x37)}
  x38<-try(ergm(g~nonzero+sum+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x38)!="try-error")){t_AICs[38]<-AIC(x38)}
  x39<-try(ergm(g~nonzero+sum+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x39)!="try-error")){t_AICs[39]<-AIC(x39)}
  x40<-try(ergm(g~nonzero+sum+nodeocov('meanRR.mm.')+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x40)!="try-error")){t_AICs[40]<-AIC(x40)}
  x41<-try(ergm(g~nonzero+sum+nodeocov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x41)!="try-error")){t_AICs[41]<-AIC(x41)}
  x42<-try(ergm(g~nonzero+sum+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x42)!="try-error")){t_AICs[42]<-AIC(x42)}
  x43<-try(ergm(g~nonzero+sum+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x43)!="try-error")){t_AICs[43]<-AIC(x43)}
  x44<-try(ergm(g~nonzero+sum+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x44)!="try-error")){t_AICs[44]<-AIC(x44)}
  x45<-try(ergm(g~nonzero+sum+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x45)!="try-error")){t_AICs[45]<-AIC(x45)}
  x46<-try(ergm(g~nonzero+sum+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x46)!="try-error")){t_AICs[46]<-AIC(x46)}
  x47<-try(ergm(g~nonzero+sum+nodeocov('rs_moy')+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x47)!="try-error")){t_AICs[47]<-AIC(x47)}
  x48<-try(ergm(g~nonzero+sum+nodeocov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x48)!="try-error")){t_AICs[48]<-AIC(x48)}
  x49<-try(ergm(g~nonzero+sum+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x49)!="try-error")){t_AICs[49]<-AIC(x49)}
  x50<-try(ergm(g~nonzero+sum+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x50)!="try-error")){t_AICs[50]<-AIC(x50)}
  x51<-try(ergm(g~nonzero+sum+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x51)!="try-error")){t_AICs[51]<-AIC(x51)}
  x52<-try(ergm(g~nonzero+sum+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x52)!="try-error")){t_AICs[52]<-AIC(x52)}
  x53<-try(ergm(g~nonzero+sum+nodeocov('MeanDensityTurtle')+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x53)!="try-error")){t_AICs[53]<-AIC(x53)}
  x54<-try(ergm(g~nonzero+sum+nodeocov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x54)!="try-error")){t_AICs[54]<-AIC(x54)}
  x55<-try(ergm(g~nonzero+sum+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x55)!="try-error")){t_AICs[55]<-AIC(x55)}
  x56<-try(ergm(g~nonzero+sum+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x56)!="try-error")){t_AICs[56]<-AIC(x56)}
  x57<-try(ergm(g~nonzero+sum+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x57)!="try-error")){t_AICs[57]<-AIC(x57)}
  x58<-try(ergm(g~nonzero+sum+nodeocov('MeanTotalFreq')+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x58)!="try-error")){t_AICs[58]<-AIC(x58)}
  x59<-try(ergm(g~nonzero+sum+nodeocov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x59)!="try-error")){t_AICs[59]<-AIC(x59)}
  x60<-try(ergm(g~nonzero+sum+nodeocov('nb_bull')+nodeocov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x60)!="try-error")){t_AICs[60]<-AIC(x60)}
  x61<-try(ergm(g~nonzero+sum+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x61)!="try-error")){t_AICs[61]<-AIC(x61)}
  x62<-try(ergm(g~nonzero+sum+nodeocov('nb_bull')+nodeocov('nb_tiger')+nodeocov('TP_Meme')+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x62)!="try-error")){t_AICs[62]<-AIC(x62)}
  x63<-try(ergm(g~nonzero+sum+nodeocov('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x63)!="try-error")){t_AICs[63]<-AIC(x63)}
  x64<-try(ergm(g~nonzero+sum+nodeocov('nb_tiger')+nodeocov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x64)!="try-error")){t_AICs[64]<-AIC(x64)}
  x65<-try(ergm(g~nonzero+sum+nodeocov('nb_tiger')+nodeocov('TP_Meme')+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x65)!="try-error")){t_AICs[65]<-AIC(x65)}
  x66<-try(ergm(g~nonzero+sum+nodeocov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x66)!="try-error")){t_AICs[66]<-AIC(x66)}
  x67<-try(ergm(g~nonzero+sum+nodeocov('TP_Meme')+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x67)!="try-error")){t_AICs[67]<-AIC(x67)}
  x68<-try(ergm(g~nonzero+sum+nodeocov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x68)!="try-error")){t_AICs[68]<-AIC(x68)}
  
  options(warn=1)
  
  b<-which.min(t_AICs)
  
  #b<-which.min(sapply(list(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, #x16,
  #                             x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30#, x31,
  #                             x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45#, x46, x47,
  #                             x48, x49, x50, x51, x52, x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63,
  #                             x64, x65, x66), AIC))
  
  networks[[i]][[7]]<-formula(paste0("x",b))
  
}

# networks[[2]][[7]]

save(networks, file="ergm_nodeocov.RData")

load("ergm_nodeocov.RData")


####
### for term 'absdiff' added to the list of the networks (8th dimension)

for (i in 1:length(networks)){
  g<-networks[[i]]
  t_AICs<-rep(NA,68)
  #options(warn=2)
  x1<-try(ergm(g~nonzero+sum+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x1)=="ergm")){t_AICs[1]<-AIC(x1)}
  x2<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x2)!="try-error")){t_AICs[2]<-AIC(x2)}
  x3<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x3)!="try-error")){t_AICs[3]<-AIC(x3)}
  x4<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+absdiff('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x4)!="try-error")){t_AICs[4]<-AIC(x4)}
  x5<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x5)!="try-error")){t_AICs[5]<-AIC(x5)}
  x6<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x6)!="try-error")){t_AICs[6]<-AIC(x6)}
  x7<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x7)!="try-error")){t_AICs[7]<-AIC(x7)}
  x8<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x8)!="try-error")){t_AICs[8]<-AIC(x8)}
  x9<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x9)!="try-error")){t_AICs[9]<-AIC(x9)}
  x10<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x10)!="try-error")){t_AICs[10]<-AIC(x10)}
  x11<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x11)!="try-error")){t_AICs[11]<-AIC(x11)}
  x12<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x12)!="try-error")){t_AICs[12]<-AIC(x12)}
  x13<-try(ergm(g~nonzero+sum+absdiff('TxoccBS')+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x13)!="try-error")){t_AICs[13]<-AIC(x13)}
  x14<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x14)!="try-error")){t_AICs[14]<-AIC(x14)}
  x15<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x15)!="try-error")){t_AICs[15]<-AIC(x15)}
  x16<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x16)!="try-error")){t_AICs[16]<-AIC(x16)}
  x17<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x17)!="try-error")){t_AICs[17]<-AIC(x17)}
  x18<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x18)!="try-error")){t_AICs[18]<-AIC(x18)}
  x19<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x19)!="try-error")){t_AICs[19]<-AIC(x19)}
  x20<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x20)!="try-error")){t_AICs[20]<-AIC(x20)}
  x21<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x21)!="try-error")){t_AICs[21]<-AIC(x21)}
  x22<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x22)!="try-error")){t_AICs[22]<-AIC(x22)}
  x23<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x23)!="try-error")){t_AICs[23]<-AIC(x23)}
  x24<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x24)!="try-error")){t_AICs[24]<-AIC(x24)}
  x25<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x25)!="try-error")){t_AICs[25]<-AIC(x25)}
  x26<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x26)!="try-error")){t_AICs[26]<-AIC(x26)}
  x27<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x27)!="try-error")){t_AICs[27]<-AIC(x27)}
  x28<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x28)!="try-error")){t_AICs[28]<-AIC(x28)}
  x29<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x29)!="try-error")){t_AICs[29]<-AIC(x29)}
  x30<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x30)!="try-error")){t_AICs[30]<-AIC(x30)}
  x31<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x31)!="try-error")){t_AICs[31]<-AIC(x31)}
  x32<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x32)!="try-error")){t_AICs[32]<-AIC(x32)}
  x33<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x33)!="try-error")){t_AICs[33]<-AIC(x33)}
  x34<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x34)!="try-error")){t_AICs[34]<-AIC(x34)}
  x35<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x35)!="try-error")){t_AICs[35]<-AIC(x35)}
  x36<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x36)!="try-error")){t_AICs[36]<-AIC(x36)}
  x37<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x37)!="try-error")){t_AICs[37]<-AIC(x37)}
  x38<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x38)!="try-error")){t_AICs[38]<-AIC(x38)}
  x39<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x39)!="try-error")){t_AICs[39]<-AIC(x39)}
  x40<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x40)!="try-error")){t_AICs[40]<-AIC(x40)}
  x41<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x41)!="try-error")){t_AICs[41]<-AIC(x41)}
  x42<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x42)!="try-error")){t_AICs[42]<-AIC(x42)}
  x43<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x43)!="try-error")){t_AICs[43]<-AIC(x43)}
  x44<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x44)!="try-error")){t_AICs[44]<-AIC(x44)}
  x45<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x45)!="try-error")){t_AICs[45]<-AIC(x45)}
  x46<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x46)!="try-error")){t_AICs[46]<-AIC(x46)}
  x47<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x47)!="try-error")){t_AICs[47]<-AIC(x47)}
  x48<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x48)!="try-error")){t_AICs[48]<-AIC(x48)}
  x49<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x49)!="try-error")){t_AICs[49]<-AIC(x49)}
  x50<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x50)!="try-error")){t_AICs[50]<-AIC(x50)}
  x51<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x51)!="try-error")){t_AICs[51]<-AIC(x51)}
  x52<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x52)!="try-error")){t_AICs[52]<-AIC(x52)}
  x53<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x53)!="try-error")){t_AICs[53]<-AIC(x53)}
  x54<-try(ergm(g~nonzero+sum+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x54)!="try-error")){t_AICs[54]<-AIC(x54)}
  x55<-try(ergm(g~nonzero+sum+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x55)!="try-error")){t_AICs[55]<-AIC(x55)}
  x56<-try(ergm(g~nonzero+sum+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x56)!="try-error")){t_AICs[56]<-AIC(x56)}
  x57<-try(ergm(g~nonzero+sum+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x57)!="try-error")){t_AICs[57]<-AIC(x57)}
  x58<-try(ergm(g~nonzero+sum+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x58)!="try-error")){t_AICs[58]<-AIC(x58)}
  x59<-try(ergm(g~nonzero+sum+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x59)!="try-error")){t_AICs[59]<-AIC(x59)}
  x60<-try(ergm(g~nonzero+sum+absdiff('nb_bull')+absdiff('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x60)!="try-error")){t_AICs[60]<-AIC(x60)}
  x61<-try(ergm(g~nonzero+sum+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x61)!="try-error")){t_AICs[61]<-AIC(x61)}
  x62<-try(ergm(g~nonzero+sum+absdiff('nb_bull')+absdiff('nb_tiger')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x62)!="try-error")){t_AICs[62]<-AIC(x62)}
  x63<-try(ergm(g~nonzero+sum+absdiff('nb_tiger')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x63)!="try-error")){t_AICs[63]<-AIC(x63)}
  x64<-try(ergm(g~nonzero+sum+absdiff('nb_tiger')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x64)!="try-error")){t_AICs[64]<-AIC(x64)}
  x65<-try(ergm(g~nonzero+sum+absdiff('nb_tiger')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x65)!="try-error")){t_AICs[65]<-AIC(x65)}
  x66<-try(ergm(g~nonzero+sum+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x66)!="try-error")){t_AICs[66]<-AIC(x66)}
  x67<-try(ergm(g~nonzero+sum+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x67)!="try-error")){t_AICs[67]<-AIC(x67)}
  x68<-try(ergm(g~nonzero+sum+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMC.compress=FALSE,MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x68)!="try-error")){t_AICs[68]<-AIC(x68)}
  
  options(warn=1)
  
  b<-which.min(t_AICs)
  
  #b<-which.min(sapply(list(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, #x16,
  #                             x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30#, x31,
  #                             x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45#, x46, x47,
  #                             x48, x49, x50, x51, x52, x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63,
  #                             x64, x65, x66), AIC))
  
  networks[[i]][[8]]<-formula(paste0("x",b))
  
}

# networks[[2]][[8]]

save(networks, file="ergm_absdiff.RData")

#load("ergm_absdiff.RData")


######################################################################################
############# STEP 5 ##################################################
#######################################################################
library(asbio)
load("ergm_absdiff.RData")
#load("ergm_absdiff.RData")
#get(load("ergm_selection1.RData"))

#get(load("ergm_selection.RData"))[[1]][[6]]

############## Extract coefficients #############

#create a table to put the results

# name | months | term | Coefficients | s.d. | p-value |
#ergm.count.tab<-as.data.frame(matrix(NA,nrow=length(networks)*11,ncol=7))
ergm.count.tab<-as.data.frame(matrix(NA,nrow=0,ncol=7))

#table to store convergence diagnostics
geweke_tab<-matrix(NA,nr=length(networks),nc=3)
rhs_tab<-matrix(NA,nr=length(networks),nc=3)

#colnames(ergm.count.tab)<-c("Shark","Month","Term","Variable","Coefficient","s.d.","p-value")
for (i in 1:length(networks)){
  
  # nodeicov
  f1<-networks[[i]][[6]]
  g<-networks[[i]]
  formula1<-paste("g~",f1[3],",response='weight', reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type='penalized')",sep="")
  md1<-eval(parse(text=paste('model=ergm(',formula1,')',sep='')))
  s1 <- summary(md1)
  #coefficient.names <- rownames(s1$coefs)
  #coefficients <- s1$coefs[, 1]
  #standard.errors <- s1$coefs[, 2]
  #significance <- s1$coefs[, 4]
  x1<-cbind(rownames(s1$coefs),s1$coefs[, 1],s1$coefs[, 2],s1$coefs[, 4])
  x11<-cbind(rep(names[[i]], nrow(x1)),rep(months[[i]], nrow(x1)),rep("nodeicov", nrow(x1)), x1)
  colnames(x11) <- colnames(ergm.count.tab)
  ergm.count.tab<-rbind(ergm.count.tab,x11)
  
  #nodeocov
  f2<-networks[[i]][[7]]
  formula2<-paste("g~",f2[3],",response='weight', reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type='penalized')",sep="")
  md2<-eval(parse(text=paste('model=ergm(',formula2,')',sep='')))
  s2 <- summary(md2)
  #coefficient.names <- rownames(s2$coefs)
  #coefficients <- s2$coefs[, 1]
  #standard.errors <- s2$coefs[, 2]
  #significance <- s2$coefs[, 4]
  x2<-cbind(rownames(s2$coefs),s2$coefs[, 1],s2$coefs[, 2],s2$coefs[, 4])
  x22<-cbind(rep(names[[i]], nrow(x2)),rep(months[[i]], nrow(x2)),rep("nodeocov", nrow(x2)), x2)
  colnames(x22) <- colnames(ergm.count.tab)
  ergm.count.tab<-rbind(ergm.count.tab,x22)
  
  #absdiff
  f3<-networks[[i]][[8]]
  formula3<-paste("g~",f3[3],",response='weight', reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type='penalized')",sep="")
  md3<-eval(parse(text=paste('model=ergm(',formula3,')',sep='')))
  s3 <- summary(md3)
  x3<-cbind(rownames(s3$coefs),s3$coefs[, 1],s3$coefs[, 2],s3$coefs[, 4])
  x33<-cbind(rep(names[[i]], nrow(x3)),rep(months[[i]], nrow(x3)),rep("absdiff", nrow(x3)), x3)
  colnames(x33) <- colnames(ergm.count.tab)
  ergm.count.tab<-rbind(ergm.count.tab,x33)
  
  mds<-list(md1,md2,md3)
  
  for(j in 1:3){
    
    mdt<-mds[[j]]
    
    #Calculate burn-in diganostics
    gt<-geweke.diag(mdt$sample,frac1=0.1,frac2=0.5)
    #our test statistic is the proportion of p values below 0.05
    geweke_tab[i,j]<-sum(pnorm(-abs(gt$z))<0.05)/length(gt$z)
    
    #Calculate Rhat diagnostics f
    #Compare first and second halves of each chain
    rhs<-numeric()
    for(ii in 1:ncol(mdt$sample)){
      Ch<-cbind(mdt$sample[1:(nrow(mdt$sample)/2),ii],mdt$sample[((nrow(mdt$sample)/2)+1):nrow(mdt$sample),ii])
      rhs[ii]<-R.hat(Ch,burn.in=0)
    }
    #our test statistic is the number of Rhat values above 1.1 (standard threshold to judge convergence) 
    rhs_tab[i,j]<-sum(rhs>1.1)
    
  }
  
}

colnames(ergm.count.tab)<-c("Shark","Month","Term","Variable","Coefficient","s.d.","p-value")
ergm.count.tab

save(ergm.count.tab, file="ergm.count.tab.RData")

#Added here to save this output 
colnames(geweke_tab)<-colnames(rhs_tab)<-c("nodeicov","nodeocov","absdiff")
save(geweke_tab, file="geweke_tab.RData")
save(rhs_tab, file="rhs_tab.RData")

#Based on the few I've been able to run it looks like the Geweke diagnostic is the most likely to be problematic. I think the best policy is to visually inspect any models where it is >0.2 - it is fairly sensitive.
#This can be done using mcmc.diagnostics(<MODELNAME>,vars.per.page=Y)
#where Y needs to be high enough to see all of them
#Here you will be able to see traces of the Markov chains and also the overall Geweke diagnostic


######################################################################################
############# STEP 6 ##################################################
#######################################################################

####### SUMMARIZE AND PLOT RESULTS ########

load("ergm.count.tab.RData")

ergm.count.tab2<-ergm.count.tab
ergm.count.tab2$Variables<-str_replace_all(ergm.count.tab2$Variable, c("nodeicov.sum.MeanHsig" = "Swell", "nodeocov.sum.MeanHsig" = "Swell", "absdiff.sum.MeanHsig" = "Swell",
                                                                       "nodeicov.sum.TxoccBS" = "Occupancy Bull Shark", "nodeocov.sum.TxoccBS" = "Occupancy Bull Shark", "absdiff.sum.TxoccBS" = "Occupancy Bull Shark",
                                                                       "nodeicov.sum.TxoccTS" = "Occupancy Tiger Shark", "nodeocov.sum.TxoccTS" = "Occupancy Tiger Shark", "absdiff.sum.TxoccTS" = "Occupancy Tiger Shark",
                                                                       "nodeicov.sum.meanRR.mm." = "Rain", "nodeocov.sum.meanRR.mm." = "Rain", "absdiff.sum.meanRR.mm." = "Rain",
                                                                       "nodeicov.sum.MeanDensityTurtle" = "Turtle density", "nodeocov.sum.MeanDensityTurtle" = "Turtle density", "absdiff.sum.MeanDensityTurtle" = "Turtle density",
                                                                       "nodeicov.sum.MeanTotalFreq" = "Human activities", "nodeocov.sum.MeanTotalFreq" = "Human activities", "absdiff.sum.MeanTotalFreq" = "Human activities",
                                                                       "nodeicov.sum.rs_moy" = "Turbidity", "nodeocov.sum.rs_moy" = "Turbidity", "absdiff.sum.rs_moy" = "Turbidity",
                                                                       "nodeicov.sum.nb_bull" = "Bull Shark abundance", "nodeocov.sum.nb_bull" = "Bull Shark abundance", "absdiff.sum.nb_bull" = "Bull Shark abundance",
                                                                       "nodeicov.sum.nb_tiger" = "Tiger Shark abundance", "nodeocov.sum.nb_tiger" = "Tiger Shark abundance", "absdiff.sum.nb_tiger" = "Tiger Shark abundance",
                                                                       "nodeicov.sum.TP_Meme" = "Same sex CRT", "nodeocov.sum.TP_Meme" = "Same sex CRT", "absdiff.sum.TP_Meme" = "Same sex CRT",
                                                                       "nodeicov.sum.TP_Diff" = "Opposite sex CRT", "nodeocov.sum.TP_Diff" = "Opposite sex CRT", "absdiff.sum.TP_Diff" = "Opposite sex CRT",
                                                                       "edgecov.sum.Geodist" = "Geodistance"))
### plot ERGM coefficients with 95%CI for each variable and each term
library(plyr)
ergm.count.tab2$Coefficient<-as.numeric(as.character(ergm.count.tab2$Coefficient))

requins<-read.table("shark.name.csv", sep=",",header=TRUE, stringsAsFactors = FALSE)
colnames(requins)<-c("Shark","Sex","Species")

dataCoeff<-merge(ergm.count.tab2,requins, by="Shark")
dataCoeff$Season<-1
dataCoeff[which(as.numeric(as.character(dataCoeff$Month))<18),11]<-"Summer"
dataCoeff[which(as.numeric(as.character(dataCoeff$Month))>17),11]<-"Winter"
ergm.count.tab2<-dataCoeff

df<-ddply(ergm.count.tab2, .(Variables,Term), summarise, mean=mean(Coefficient), se = sd(Coefficient)/sqrt(length(Coefficient)),lowerCI = mean(Coefficient)-qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), upperCI = mean(Coefficient)+qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), n = length(Coefficient))


pd <- position_dodge(0.8) # move them .05 to the left and right

#ggplot(df, aes(x=factor(Variables, level = c('Swell', 'Turbidity', 'Rain','Human activities','Turtle density','Same sex CRT','Opposite sex CRT','Bull Shark abundance', 'Occupancy Bull Shark', 'Tiger Shark abundance','Occupancy Tiger Shark','Geodistance','sum', 'nonzero')), y=mean, color=Term)) + 
ggplot(df, aes(x=factor(Variables, level = c('Geodistance','sum','nonzero','Occupancy Tiger Shark','Tiger Shark abundance','Occupancy Bull Shark','Bull Shark abundance','Opposite sex CRT','Same sex CRT','Human activities','Turtle density','Rain','Turbidity','Swell')), y=mean, color=Term)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0, size=1, position=pd) +
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0, size=0.5, position=pd) +
  geom_point(aes(size=n), alpha=0.5, position=pd)+
  scale_size(breaks=c(1,10,20,30,40), range = c(1,4)) +
  geom_hline(yintercept = 0,linetype = "dotted")+
  scale_color_brewer(palette="Dark2")+
  coord_flip()+
  xlab('Variable')

ergm.count.tab2[which(ergm.count.tab2$Variables=="Opposite sex CRT"),1:5]

range(ergm.count.tab2[which(ergm.count.tab2$Variables=="Opposite sex CRT"),5])



#### Graph for Sex

df_sex<-ddply(ergm.count.tab2, .(Variables,Term,Sex), summarise, mean=mean(Coefficient), se = sd(Coefficient)/sqrt(length(Coefficient)),lowerCI = mean(Coefficient)-qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), upperCI = mean(Coefficient)+qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), n = length(Coefficient))


pd <- position_dodge(0.8) # move them .05 to the left and right
ggplot(df_sex, aes(x=factor(Variables, level = c('Geodistance','sum','nonzero','Occupancy Tiger Shark','Tiger Shark abundance','Occupancy Bull Shark','Bull Shark abundance','Opposite sex CRT','Same sex CRT','Human activities','Turtle density','Rain','Turbidity','Swell')), y=mean, color=Term)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0, size=1, position=pd) +
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0, size=0.5, position=pd) +
  geom_point(aes(size=n), alpha=0.5, position=pd)+
  scale_size(breaks=c(1,10,20,30,40), range = c(1,4)) +
  geom_hline(yintercept = 0,linetype = "dotted")+
  scale_color_brewer(palette="Dark2")+
  coord_flip()+
  facet_grid(~Sex)+
  xlab('Variable')


# Graph for Season

df_season<-ddply(ergm.count.tab2, .(Variables, Term, Season), summarise, mean=mean(Coefficient), se = sd(Coefficient)/sqrt(length(Coefficient)),lowerCI = mean(Coefficient)-qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), upperCI = mean(Coefficient)+qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), n = length(Coefficient))

pd <- position_dodge(0.8) # move them .05 to the left and right
ggplot(df_season, aes(x=factor(Variables, level = c('Geodistance','sum','nonzero','Occupancy Tiger Shark','Tiger Shark abundance','Occupancy Bull Shark','Bull Shark abundance','Opposite sex CRT','Same sex CRT','Human activities','Turtle density','Rain','Turbidity','Swell')), y=mean, color=Term)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0, size=1, position=pd) +
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0, size=0.5, position=pd) +
  geom_point(aes(size=n), alpha=0.5, position=pd)+
  scale_size(breaks=c(1,10,20,30,40), range = c(1,4)) +
  geom_hline(yintercept = 0,linetype = "dotted")+
  scale_color_brewer(palette="Dark2")+
  coord_flip()+
  facet_grid(~Season)+
  xlab('Variable')







########## Plot Movement network on a map

#load networks
#load("malist.RData")
load("network_list.RData")
load("ergm_selection.RData")

# Sum number of mouvements for all months & all inividuals
rows <- dim(mouvments)[1]
cols <- dim(mouvments)[2]

result <- matrix(0, nrow=rows, ncol=cols)
for (i in seq(rows)) {
  for (j in seq(cols)) {
    result[i,j] <- sum(mouvments[i,j,,])
  }
}

rownames(result)<-rownames(mouvments)
colnames(result)<-colnames(mouvments)
aggregated<-result


#agg<-as.network(aggregated, matrix.type = "adjacency",directed = FALSE,names.eval = "weight",ignore.eval = FALSE,loops=F)

#plot(agg)

aggregated2<-aggregated
aggregated2[aggregated2<15]<-0     ### only edges with more than 10 movements are indicated 
colnames(aggregated2)<-colnames(mouvments)
rownames(aggregated2)<-rownames(mouvments)


# Coordinates of sations and associated color codes 
sta<-as.data.frame(unique(tableIbis$code_lieu))
colnames(sta)<-"code_lieu"

sta2<-unique(bdd[,c("lat","lon","code_lieu")])
sta2<-sta2[order(sta2$code_lieu),]
sta2<-sta2[-6,]

XY_arrays<-merge(sta,sta2,by.x="code_lieu", all.x=TRUE)
#XY_arrays<-unique(bdd[,c("lat","lon","code_lieu")])
XY_arrays<-unique(XY_arrays[which(XY_arrays[,1] %in% station.names),])
rownames(XY_arrays)<-XY_arrays[,1]
XY_arrays[15,]<-c("CM-LSGS",-21.06685,55.1936)
XY_arrays[,2]<-as.numeric(XY_arrays[,2])
XY_arrays[,3]<-as.numeric(XY_arrays[,3])

#z["STROSE"]<-c(-21.36)

island <- getData('GADM', country='REU', level=0)
#plot(island,col="gray",main="Carte de la Reunion",cex.main=2)

diag(aggregated2)<-0
tot<-sum(apply(aggregated2, 2, sum))
aggregated3<-aggregated2/tot

###### Function to create vector of variable to color
val2col<-function(z, zlim, col = chlPal(10), breaks){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break
                                               than colour")}
    }
  if(missing(breaks) & !missing(zlim)){
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  CUT <- cut(z, breaks=breaks,include.lowest=T)
  colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
  return(colorlevels)
  }

#####


net<-as.network(aggregated2, matrix.type = "adjacency", directed = FALSE,names.eval = "weight",ignore.eval = FALSE,loops=F)
net %v% "vertex.names"<-rownames(aggregated2)
net %e% "weight"<-(net %e% "weight")

nc = 60  #number of colors
pal <- wes_palette("Zissou1", nc, type = "continuous")
z=net %e% "weight"
pal2=val2col(z=z,zlim=max(z),pal,breaks=seq(0,max(z),length.out = nc+1))
rownames(XY_arrays)=XY_arrays[,1]
par(bg="gray35")
plot(island,col="gray25")
plot.network(net,coord=XY_arrays[net %v% 'vertex.names',3:2],vertex.col="white",vertex.cex=0.3,vertex.size=0.3, arrowhead.cex = 0.3, edge.lw=(net %e% "weight")/100, edge.col=pal2,  usecurve = FALSE, edge.curve = 1e-03,uselen= FALSE,edge.len = 0.01,suppress.axes=FALSE,new=FALSE)
addscalebar(linecol = "white",label.col="white")+
  addnortharrow(border="white",text.col = "white", scale = 0.7)+
  gradientLegend(valRange=range(z),pos=0.5, side=4,color=pal)


#dev.off()


write.csv(aggregated2,"movements.csv")

####
sub<-unique(ergm.count.tab[,1:2])

for i in unique(sub$Shark) & j in unique(sub$Shark)
list(mouvments[,,unique(ergm.count.tab[,1:2])])

result <- matrix(0, nrow=rows, ncol=cols)
for (i in unique(sub$Shark)) {
  for (j in unique(sub$Month)) {
    result[i,j] <- sum(mouvments[i,j,,])
  }
}


matsub<-list()
for (i in unique(sub$Shark)) {
  for (j in unique(sub$Month))
  matsub<-mouvments[,,i,j]
}
mouvments[,,"Adelante","1"]

