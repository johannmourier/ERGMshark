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

# Set the working directory 
setwd("/Users/johannmourier/Documents/Scripts/CHARC/Stage_Angelique/Etape 7 - Modele ERGM/NEW")


######################################################################################
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

#save(param1,file="Param scale3.RData")

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



#############################################################################
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
          if(network.edgecount(g) >=5){
            
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

##################################################################
####################  RUN ERGM Models  ###########################

#table storing ergm coefficients in a table
ergm.count.tab<-as.data.frame(matrix(0,nrow=length(malist),ncol=77))
for (i in 1:length(malist)){ 
  
    ergm.count.tab[i,1]<-names[[i]]
    ergm.count.tab[i,2]<-months[[i]]
    g<-networks[[i]]
    x<-ergm(g~nonzero+nodecov('TxoccBS')+diff('TxoccBS')+nodecov('TxoccTS')+diff('TxoccTS')+nodecov('MeanHsig')+diff('MeanHsig')+nodecov('MeanDensityTurtle')+diff('MeanDensityTurtle')+nodecov('meanRR.mm.')+diff('meanRR.mm.')+nodecov('MeanTotalFreq')+diff('MeanTotalFreq')+nodecov('rs_moy')+diff('rs_moy')+nodecov('nb_bull')+diff('nb_bull')+nodecov('nb_tiger')+diff('nb_tiger')+nodecov('TP_Meme')+diff('TP_Meme')+nodecov('TP_Diff')+diff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized"))
    #summary(x)
    ergm.count.tab[i,3:5]<-(summary(x)$coefs)[1,-3]    #nonzero
    ergm.count.tab[i,6:8]<-(summary(x)$coefs)[2,-3]    #nodecov.sum.TxoccBS
    ergm.count.tab[i,9:11]<-(summary(x)$coefs)[3,-3]   #diff.sum.t-h.TxoccBS
    ergm.count.tab[i,12:14]<-(summary(x)$coefs)[4,-3]  #nodecov.sum.TxoccTS
    ergm.count.tab[i,15:17]<-(summary(x)$coefs)[5,-3]  #diff.sum.t-h.TxoccTS
    ergm.count.tab[i,18:20]<-(summary(x)$coefs)[6,-3]  #nodecov.sum.MeanHsig
    ergm.count.tab[i,21:23]<-(summary(x)$coefs)[7,-3]  #diff.sum.t-h.MeanHsig
    ergm.count.tab[i,24:26]<-(summary(x)$coefs)[8,-3]  #nodecov.sum.MeanDensityTurtle
    ergm.count.tab[i,27:29]<-(summary(x)$coefs)[9,-3]  #diff.sum.t-h.MeanDensityTurtle
    ergm.count.tab[i,30:32]<-(summary(x)$coefs)[10,-3] #nodecov.sum.meanRR.mm.
    ergm.count.tab[i,33:35]<-(summary(x)$coefs)[11,-3] #diff.sum.t-h.meanRR.mm.
    ergm.count.tab[i,36:38]<-(summary(x)$coefs)[12,-3] #nodecov.sum.MeanTotalFreq
    ergm.count.tab[i,39:41]<-(summary(x)$coefs)[13,-3] #diff.sum.t-h.MeanTotalFreq
    ergm.count.tab[i,42:44]<-(summary(x)$coefs)[14,-3] #nodecov.sum.rs_moy
    ergm.count.tab[i,45:47]<-(summary(x)$coefs)[15,-3] #diff.sum.t-h.rs_moy
    ergm.count.tab[i,48:50]<-(summary(x)$coefs)[16,-3] #nodecov.sum.nb_bull
    ergm.count.tab[i,51:53]<-(summary(x)$coefs)[17,-3] #diff.sum.t-h.nb_bull
    ergm.count.tab[i,54:56]<-(summary(x)$coefs)[18,-3] #nodecov.sum.nb_tiger
    ergm.count.tab[i,57:59]<-(summary(x)$coefs)[19,-3] #diff.sum.t-h.nb_tiger
    ergm.count.tab[i,60:62]<-(summary(x)$coefs)[20,-3] #nodecov.sum.TP_Meme
    ergm.count.tab[i,63:65]<-(summary(x)$coefs)[21,-3] #diff.sum.t-h.TP_Meme
    ergm.count.tab[i,66:68]<-(summary(x)$coefs)[22,-3] #nodecov.sum.TP_Diff
    ergm.count.tab[i,69:71]<-(summary(x)$coefs)[23,-3] #diff.sum.t-h.TP_Diff
    ergm.count.tab[i,72:74]<-(summary(x)$coefs)[24,-3] #edgecov.sum.Geodist
    ergm.count.tab[i,75]<-network.edgecount(g)
    ergm.count.tab[i,76]<-length(network.vertex.names(g))
    ergm.count.tab[i,77]<-sum((g %e% 'weight')==1)
}

colnames(ergm.count.tab)<- c("Requin","Mois","nonzero_Estimate","nonzero_SE","nonzero_pv","TxoccBS_Estimate","TxoccBS_SE","TxoccBS_pv",
                             "TxoccTS_Estimate","TxoccTS_SE","TxoccTS_pv","meanRR.mm._Estimate","meanRR.mm._SE","meanRR.mm._pv",
                             "MeanHsig_Estimate","MeanHsig_SE","MeanHsig_pv",
                             "MeanDensityTurtle_Estimate","MeanDensityTurtle_SE","MeanDensityTurtle_pv","MeanTotalFreq_Estimate","MeanTotalFreq_SE","MeanTotalFreq_pv",
                             "rs_moy_Estimate","rs_moy_SE","rs_moy_pv",
                             "nb_bull_Estimate","nb_bull_SE","nb_bull_pv","nb_tiger_Estimate","nb_tiger_SE","nb_tiger_pv",
                             "TP_Meme_Estimate","TP_Meme_SE","TP_Meme_pv","TP_Diff_Estimate","TP_Diff_SE","TP_Diff_pv",
                             "Geodist_Estimate","Geodist_SE","Geodist_pv","Nb de edges","Nb de noeuds","weightEqual1")

ergm.count.tab


