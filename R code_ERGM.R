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


#### create matrix : months / stations / sharks / residency with same-different sex
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

#########################     NETWORK SELECTION      ##########################

############################
#############
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
#This is to create the order in the original step 3
node_atts<-c("TxoccBS","TxoccTS","meanRR.mm.","MeanHsig","MeanDensityTurtle",
             "MeanTotalFreq","rs_moy","nb_bull","nb_tiger","TP_Meme","TP_Diff")
####### Procedure that selects the terms of the ergm model for each shark and each month
for (i in shark.names ){
  date<-colnames(requins)[which(requins[i,]>time.cutoff)] # month in which the shark is tagged at least 20 days
  if (length(date)>0){
    for (j in date){
      stat<-rownames(stations)[which(as.numeric(as.character(stations[,j]))>time.cutoff)] #stations deployed more than 20 days during month j
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
          for(ii in 1:length(node_atts)){
            #colSums here is an arbitrary choice, could equally easily be rowSums. Or
            #we could use both....
            t_sds<-aggregate(ver.attr1[,1],by=list(sign(colSums(as.matrix(g))+rowSums(as.matrix(g)))),sd)[,2]
            sum(is.na(t_sds))
            
            if(sum(is.na(t_sds))==0&sum(t_sds==0)==0){
              print(paste("networks[[",l,"]] %v% '",node_atts[ii],"' <- ver.attr1[,",ii,"]",sep=''))
              eval(parse(text=paste("g %v% '",node_atts[ii],"' <- ver.attr1[,",ii,"]",sep='')))
            }
          }
          g %e% "Geodist"<-Geodist[rownames(ver.attr1),rownames(ver.attr1)] 
          ## we keep more than 5 edges
          if(network.edgecount(g) >=10){
            
            networks[[l]]<-g 
            names[[l]]<-i      
            months[[l]]<-j      
            mylist <- list( networks , names , months)
            l<- l+1
            rm(list = c("stat","env","g","ver.attr","ver.attr1"))
            
            print(paste("Le requin",i,"mois",j))
            
          }
        }
        
      }
    }
  }
  rm(list="date")
}

#Have a quick check of months
table(unlist(months))
sum(table(unlist(months)))


######################################################################################
############# STEP 3 ##################################################
#######################################################################

#################### MODEL SELECTION ########################
#############################################################

### for each network we select the best model from smallest AIC value

####
### for term 'nodecov' added to the list of the networks (6th dimension)
for (i in 1:44){
  g<-networks[[i]]
  t_AICs<-rep(NA,68)
  #options(warn=2)
  x1<-try(ergm(g~nonzero+sum+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x1)=="ergm")){
    t_ms<-summary(x1)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[1]<-AIC(x1)
    }
  }
  x2<-try(ergm(g~nonzero+sum+nodecov('TxoccTS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.steplength.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x2)=="ergm")){
    t_ms<-summary(x2)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[2]<-AIC(x2)
    }
  }
  x3<-try(ergm(g~nonzero+sum+nodecov('TxoccTS')+nodecov('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x3)=="ergm")){
    t_ms<-summary(x3)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[3]<-AIC(x3)
    }
  }
  x4<-try(ergm(g~nonzero+sum+nodecov('TxoccTS')+nodecov('MeanHsig')+nodecov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x4)=="ergm")){
    t_ms<-summary(x4)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[4]<-AIC(x4)
    }
  }
  x5<-try(ergm(g~nonzero+sum+nodecov('TxoccTS')+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x5)=="ergm")){
    t_ms<-summary(x5)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[5]<-AIC(x5)
    }
  }
  x6<-try(ergm(g~nonzero+sum+nodecov('TxoccTS')+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x6)=="ergm")){
    t_ms<-summary(x6)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[6]<-AIC(x6)
    }
  }
  x7<-try(ergm(g~nonzero+sum+nodecov('TxoccTS')+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x7)=="ergm")){
    t_ms<-summary(x7)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[7]<-AIC(x7)
    }
  }
  x8<-try(ergm(g~nonzero+sum+nodecov('TxoccTS')+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x8)=="ergm")){
    t_ms<-summary(x8)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[8]<-AIC(x8)
    }
  }
  x9<-try(ergm(g~nonzero+sum+nodecov('TxoccTS')+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x9)=="ergm")){
    t_ms<-summary(x9)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[9]<-AIC(x9)
    }
  }
  x10<-try(ergm(g~nonzero+sum+nodecov('TxoccTS')+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x10)=="ergm")){
    t_ms<-summary(x10)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[10]<-AIC(x10)
    }
  }
  x11<-try(ergm(g~nonzero+sum+nodecov('TxoccTS')+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+nodecov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x11)=="ergm")){
    t_ms<-summary(x11)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[11]<-AIC(x11)
    }
  }
  x12<-try(ergm(g~nonzero+sum+nodecov('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x12)=="ergm")){
    t_ms<-summary(x12)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[12]<-AIC(x12)
    }
  }
  x13<-try(ergm(g~nonzero+sum+nodecov('MeanHsig')+nodecov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x13)=="ergm")){
    t_ms<-summary(x13)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[13]<-AIC(x13)
    }
  }
  x14<-try(ergm(g~nonzero+sum+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x14)=="ergm")){
    t_ms<-summary(x14)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[14]<-AIC(x14)
    }
  }
  x15<-try(ergm(g~nonzero+sum+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x15)=="ergm")){
    t_ms<-summary(x15)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[15]<-AIC(x15)
    }
  }
  x16<-try(ergm(g~nonzero+sum+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x16)=="ergm")){
    t_ms<-summary(x16)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[16]<-AIC(x16)
    }
  }
  x17<-try(ergm(g~nonzero+sum+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x17)=="ergm")){
    t_ms<-summary(x17)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[17]<-AIC(x17)
    }
  }
  x18<-try(ergm(g~nonzero+sum+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x18)=="ergm")){
    t_ms<-summary(x18)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[18]<-AIC(x18)
    }
  }
  x19<-try(ergm(g~nonzero+sum+nodecov('MeanHsig')+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+nodecov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x19)=="ergm")){
    t_ms<-summary(x19)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[19]<-AIC(x19)
    }
  }
  x20<-try(ergm(g~nonzero+sum+nodecov('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x20)=="ergm")){
    t_ms<-summary(x20)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[20]<-AIC(x20)
    }
  }
  x21<-try(ergm(g~nonzero+sum+nodecov('meanRR.mm.')+nodecov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x21)=="ergm")){
    t_ms<-summary(x21)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[21]<-AIC(x21)
    }
  }
  x22<-try(ergm(g~nonzero+sum+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x22)=="ergm")){
    t_ms<-summary(x22)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[22]<-AIC(x22)
    }
  }
  x23<-try(ergm(g~nonzero+sum+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x23)=="ergm")){
    t_ms<-summary(x23)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[23]<-AIC(x23)
    }
  }
  x24<-try(ergm(g~nonzero+sum+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x24)=="ergm")){
    t_ms<-summary(x24)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[24]<-AIC(x24)
    }
  }
  x25<-try(ergm(g~nonzero+sum+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x25)=="ergm")){
    t_ms<-summary(x25)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[25]<-AIC(x25)
    }
  }
  x26<-try(ergm(g~nonzero+sum+nodecov('meanRR.mm.')+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+nodecov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x26)=="ergm")){
    t_ms<-summary(x26)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[26]<-AIC(x26)
    }
  }
  x27<-try(ergm(g~nonzero+sum+nodecov('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x27)=="ergm")){
    t_ms<-summary(x27)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[27]<-AIC(x27)
    }
  }
  x28<-try(ergm(g~nonzero+sum+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x28)=="ergm")){
    t_ms<-summary(x28)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[28]<-AIC(x28)
    }
  }
  x29<-try(ergm(g~nonzero+sum+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x29)=="ergm")){
    t_ms<-summary(x29)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[29]<-AIC(x29)
    }
  }
  x30<-try(ergm(g~nonzero+sum+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x30)=="ergm")){
    t_ms<-summary(x30)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[30]<-AIC(x30)
    }
  }
  x31<-try(ergm(g~nonzero+sum+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x31)=="ergm")){
    t_ms<-summary(x31)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[31]<-AIC(x31)
    }
  }
  x32<-try(ergm(g~nonzero+sum+nodecov('rs_moy')+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+nodecov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x32)=="ergm")){
    t_ms<-summary(x32)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[32]<-AIC(x32)
    }
  }
  x33<-try(ergm(g~nonzero+sum+nodecov('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x33)=="ergm")){
    t_ms<-summary(x33)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[33]<-AIC(x33)
    }
  }
  x34<-try(ergm(g~nonzero+sum+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x34)=="ergm")){
    t_ms<-summary(x34)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[34]<-AIC(x34)
    }
  }
  x35<-try(ergm(g~nonzero+sum+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x35)=="ergm")){
    t_ms<-summary(x35)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[35]<-AIC(x35)
    }
  }
  x36<-try(ergm(g~nonzero+sum+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x36)=="ergm")){
    t_ms<-summary(x36)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[36]<-AIC(x36)
    }
  }
  x37<-try(ergm(g~nonzero+sum+nodecov('MeanDensityTurtle')+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+nodecov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x37)=="ergm")){
    t_ms<-summary(x37)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[37]<-AIC(x37)
    }
  }
  x38<-try(ergm(g~nonzero+sum+nodecov('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x38)=="ergm")){
    t_ms<-summary(x38)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[38]<-AIC(x38)
    }
  }
  x39<-try(ergm(g~nonzero+sum+nodecov('MeanTotalFreq')+nodecov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x39)=="ergm")){
    t_ms<-summary(x39)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[39]<-AIC(x39)
    }
  }
  x40<-try(ergm(g~nonzero+sum+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x40)=="ergm")){
    t_ms<-summary(x40)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[40]<-AIC(x40)
    }
  }
  x41<-try(ergm(g~nonzero+sum+nodecov('MeanTotalFreq')+nodecov('nb_bull')+nodecov('TP_Meme')+nodecov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x41)=="ergm")){
    t_ms<-summary(x41)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[41]<-AIC(x41)
    }
  }
  x42<-try(ergm(g~nonzero+sum+nodecov('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x42)=="ergm")){
    t_ms<-summary(x42)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[42]<-AIC(x42)
    }
  }
  x43<-try(ergm(g~nonzero+sum+nodecov('nb_bull')+nodecov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x43)=="ergm")){
    t_ms<-summary(x43)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[43]<-AIC(x43)
    }
  }
  x44<-try(ergm(g~nonzero+sum+nodecov('nb_bull')+nodecov('TP_Meme')+nodecov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x44)=="ergm")){
    t_ms<-summary(x44)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[44]<-AIC(x44)
    }
  }
  x45<-try(ergm(g~nonzero+sum+nodecov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x45)=="ergm")){
    t_ms<-summary(x45)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[45]<-AIC(x45)
    }
  }
  x46<-try(ergm(g~nonzero+sum+nodecov('TP_Meme')+nodecov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x46)=="ergm")){
    t_ms<-summary(x46)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[46]<-AIC(x46)
    }
  }
  x47<-try(ergm(g~nonzero+sum+nodecov('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x47)=="ergm")){
    t_ms<-summary(x47)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[47]<-AIC(x47)
    }
  }
  x48<-try(ergm(g~nonzero+sum+nodecov('TP_Meme')+nodecov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x48)=="ergm")){
    t_ms<-summary(x48)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[48]<-AIC(x48)
    }
  }
  x49<-try(ergm(g~nonzero+sum+nodecov('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x49)=="ergm")){
    t_ms<-summary(x49)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[49]<-AIC(x49)
    }
  }
  
  options(warn=1)
  
  b<-which.min(t_AICs)
  networks[[i]][[6]]<-eval(parse(text=paste("formula(x",b,")",sep='')))
  
}


####
### for term 'absdiff' added to the list of the networks (8th dimension)
for (i in 1:44){
  g<-networks[[i]]
  t_AICs<-rep(NA,68)
  #options(warn=2)
  x1<-try(ergm(g~nonzero+sum+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x1)=="ergm")){
    t_ms<-summary(x1)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[1]<-AIC(x1)
    }
  }
  x2<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.steplength.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x2)=="ergm")){
    t_ms<-summary(x2)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[2]<-AIC(x2)
    }
  }
  x3<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x3)=="ergm")){
    t_ms<-summary(x3)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[3]<-AIC(x3)
    }
  }
  x4<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x4)=="ergm")){
    t_ms<-summary(x4)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[4]<-AIC(x4)
    }
  }
  x5<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x5)=="ergm")){
    t_ms<-summary(x5)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[5]<-AIC(x5)
    }
  }
  x6<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x6)=="ergm")){
    t_ms<-summary(x6)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[6]<-AIC(x6)
    }
  }
  x7<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x7)=="ergm")){
    t_ms<-summary(x7)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[7]<-AIC(x7)
    }
  }
  x8<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x8)=="ergm")){
    t_ms<-summary(x8)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[8]<-AIC(x8)
    }
  }
  x9<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x9)=="ergm")){
    t_ms<-summary(x9)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[9]<-AIC(x9)
    }
  }
  x10<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x10)=="ergm")){
    t_ms<-summary(x10)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[10]<-AIC(x10)
    }
  }
  x11<-try(ergm(g~nonzero+sum+absdiff('TxoccTS')+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x11)=="ergm")){
    t_ms<-summary(x11)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[11]<-AIC(x11)
    }
  }
  x12<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x12)=="ergm")){
    t_ms<-summary(x12)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[12]<-AIC(x12)
    }
  }
  x13<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x13)=="ergm")){
    t_ms<-summary(x13)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[13]<-AIC(x13)
    }
  }
  x14<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x14)=="ergm")){
    t_ms<-summary(x14)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[14]<-AIC(x14)
    }
  }
  x15<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x15)=="ergm")){
    t_ms<-summary(x15)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[15]<-AIC(x15)
    }
  }
  x16<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x16)=="ergm")){
    t_ms<-summary(x16)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[16]<-AIC(x16)
    }
  }
  x17<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x17)=="ergm")){
    t_ms<-summary(x17)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[17]<-AIC(x17)
    }
  }
  x18<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x18)=="ergm")){
    t_ms<-summary(x18)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[18]<-AIC(x18)
    }
  }
  x19<-try(ergm(g~nonzero+sum+absdiff('MeanHsig')+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x19)=="ergm")){
    t_ms<-summary(x19)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[19]<-AIC(x19)
    }
  }
  x20<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x20)=="ergm")){
    t_ms<-summary(x20)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[20]<-AIC(x20)
    }
  }
  x21<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x21)=="ergm")){
    t_ms<-summary(x21)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[21]<-AIC(x21)
    }
  }
  x22<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x22)=="ergm")){
    t_ms<-summary(x22)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[22]<-AIC(x22)
    }
  }
  x23<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x23)=="ergm")){
    t_ms<-summary(x23)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[23]<-AIC(x23)
    }
  }
  x24<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x24)=="ergm")){
    t_ms<-summary(x24)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[24]<-AIC(x24)
    }
  }
  x25<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x25)=="ergm")){
    t_ms<-summary(x25)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[25]<-AIC(x25)
    }
  }
  x26<-try(ergm(g~nonzero+sum+absdiff('meanRR.mm.')+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x26)=="ergm")){
    t_ms<-summary(x26)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[26]<-AIC(x26)
    }
  }
  x27<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x27)=="ergm")){
    t_ms<-summary(x27)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[27]<-AIC(x27)
    }
  }
  x28<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x28)=="ergm")){
    t_ms<-summary(x28)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[28]<-AIC(x28)
    }
  }
  x29<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x29)=="ergm")){
    t_ms<-summary(x29)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[29]<-AIC(x29)
    }
  }
  x30<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x30)=="ergm")){
    t_ms<-summary(x30)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[30]<-AIC(x30)
    }
  }
  x31<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x31)=="ergm")){
    t_ms<-summary(x31)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[31]<-AIC(x31)
    }
  }
  x32<-try(ergm(g~nonzero+sum+absdiff('rs_moy')+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x32)=="ergm")){
    t_ms<-summary(x32)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[32]<-AIC(x32)
    }
  }
  x33<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x33)=="ergm")){
    t_ms<-summary(x33)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[33]<-AIC(x33)
    }
  }
  x34<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x34)=="ergm")){
    t_ms<-summary(x34)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[34]<-AIC(x34)
    }
  }
  x35<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x35)=="ergm")){
    t_ms<-summary(x35)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[35]<-AIC(x35)
    }
  }
  x36<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x36)=="ergm")){
    t_ms<-summary(x36)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[36]<-AIC(x36)
    }
  }
  x37<-try(ergm(g~nonzero+sum+absdiff('MeanDensityTurtle')+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x37)=="ergm")){
    t_ms<-summary(x37)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[37]<-AIC(x37)
    }
  }
  x38<-try(ergm(g~nonzero+sum+absdiff('MeanTotalFreq')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x38)=="ergm")){
    t_ms<-summary(x38)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[38]<-AIC(x38)
    }
  }
  x39<-try(ergm(g~nonzero+sum+absdiff('MeanTotalFreq')+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x39)=="ergm")){
    t_ms<-summary(x39)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[39]<-AIC(x39)
    }
  }
  x40<-try(ergm(g~nonzero+sum+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x40)=="ergm")){
    t_ms<-summary(x40)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[40]<-AIC(x40)
    }
  }
  x41<-try(ergm(g~nonzero+sum+absdiff('MeanTotalFreq')+absdiff('nb_bull')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x41)=="ergm")){
    t_ms<-summary(x41)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[41]<-AIC(x41)
    }
  }
  x42<-try(ergm(g~nonzero+sum+absdiff('nb_bull')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x42)=="ergm")){
    t_ms<-summary(x42)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[42]<-AIC(x42)
    }
  }
  x43<-try(ergm(g~nonzero+sum+absdiff('nb_bull')+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x43)=="ergm")){
    t_ms<-summary(x43)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[43]<-AIC(x43)
    }
  }
  x44<-try(ergm(g~nonzero+sum+absdiff('nb_bull')+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x44)=="ergm")){
    t_ms<-summary(x44)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[44]<-AIC(x44)
    }
  }
  x45<-try(ergm(g~nonzero+sum+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x45)=="ergm")){
    t_ms<-summary(x45)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[45]<-AIC(x45)
    }
  }
  x46<-try(ergm(g~nonzero+sum+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x46)=="ergm")){
    t_ms<-summary(x46)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[46]<-AIC(x46)
    }
  }
  x47<-try(ergm(g~nonzero+sum+absdiff('TP_Meme')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x47)=="ergm")){
    t_ms<-summary(x47)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[47]<-AIC(x47)
    }
  }
  x48<-try(ergm(g~nonzero+sum+absdiff('TP_Meme')+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x48)=="ergm")){
    t_ms<-summary(x48)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[48]<-AIC(x48)
    }
  }
  x49<-try(ergm(g~nonzero+sum+absdiff('TP_Diff')+edgecov(Geodist,'Geodist'),response='weight',reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type="penalized")))
  if(isTRUE(class(x49)=="ergm")){
    t_ms<-summary(x49)
    if(sum(is.na(t_ms$coefs[,2]))==0){
      t_AICs[49]<-AIC(x49)
    }
  }
  
  options(warn=1)
  
  b<-which.min(t_AICs)
  networks[[i]][[7]]<-eval(parse(text=paste("formula(x",b,")",sep='')))
  
}

save(networks, file="ergm_selected_formula.RData")

#########


######################################################################################
############# STEP 4 ##################################################
#######################################################################
library(asbio)
load("ergm_selected_formula.RData")


############## Extract coefficients #############

#create a table to put the results

# name | months | term | Coefficients | s.d. | p-value |
ergm.count.tab<-as.data.frame(matrix(NA,nrow=0,ncol=7))

#table to store convergence diagnostics
geweke_tab<-matrix(NA,nr=length(networks),nc=2)
rhs_tab<-matrix(NA,nr=length(networks),nc=2)

for (i in 1:length(networks)){
  
  # nodecov
  f1<-networks[[i]][[6]]
  g<-networks[[i]]
  formula1<-paste("g~",f1[3],",response='weight', reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type='penalized')",sep="")
  md1<-eval(parse(text=paste('model=ergm(',formula1,')',sep='')))
  s1 <- summary(md1)
  x1<-cbind(rownames(s1$coefs),s1$coefs[, 1],s1$coefs[, 2],s1$coefs[, 4])
  x11<-cbind(rep(names[[i]], nrow(x1)),rep(months[[i]], nrow(x1)),rep("nodecov", nrow(x1)), x1)
  colnames(x11) <- colnames(ergm.count.tab)
  ergm.count.tab<-rbind(ergm.count.tab,x11)
  
  #absdiff
  f2<-networks[[i]][[7]]
  formula2<-paste("g~",f2[3],",response='weight', reference = ~Poisson,control=control.ergm(MCMLE.Hummel.maxit=1000,MPLE.type='penalized')",sep="")
  md2<-eval(parse(text=paste('model=ergm(',formula2,')',sep='')))
  s2 <- summary(md2)
  x2<-cbind(rownames(s2$coefs),s2$coefs[, 1],s2$coefs[, 2],s2$coefs[, 4])
  x22<-cbind(rep(names[[i]], nrow(x2)),rep(months[[i]], nrow(x2)),rep("absdiff", nrow(x2)), x2)
  colnames(x22) <- colnames(ergm.count.tab)
  ergm.count.tab<-rbind(ergm.count.tab,x22)
  
  
  mds<-list(md1,md2)
  
  for(j in 1:2){
    
    mdt<-mds[[j]]
    
    #Calculate burn-in diganostics
    gt<-geweke.diag(mdt$sample,frac1=0.1,frac2=0.5)
    #our test statistic is the proportion of p values below 0.05
    geweke_tab[i,j]<-sum(pnorm(-abs(gt[[1]]$z))*2<0.05)/length(gt[[1]]$z)
    #Calculate Rhat diagnostics f
    #Compare first and second halves of each chain
    rhs<-numeric()
    for(ii in 1:ncol(mdt$sample[[1]])){
      Ch<-cbind(mdt$sample[[1]][1:(nrow(mdt$sample[[1]])/2),ii],mdt$sample[[1]][((nrow(mdt$sample[[1]])/2)+1):nrow(mdt$sample[[1]]),ii])
      rhs[ii]<-R.hat(Ch,burn.in=0)
    }
    #our test statistic is the number of Rhat values above 1.1 (standard threshold to judge convergence) 
    rhs_tab[i,j]<-sum(rhs>1.1)
    
  }
  
  print(i)
}

colnames(ergm.count.tab)<-c("Shark","Month","Term","Variable","Coefficient","s.d.","p-value")
ergm.count.tab

save(ergm.count.tab, file="ergm.count.tab.RData")

#Added here to save this output 
colnames(geweke_tab)<-colnames(rhs_tab)<-c("nodecov","absdiff")
save(geweke_tab, file="geweke_tab.RData")
save(rhs_tab, file="rhs_tab.RData")

load("geweke_tab.RData")



######################################################################################
############# STEP 5 ##################################################
#######################################################################

####### SUMMARIZE AND PLOT RESULTS ########

load("ergm.count.tab.RData")

ergm.count.tab2<-ergm.count.tab


ergm.count.tab2$Variables<-str_replace_all(ergm.count.tab2$Variable, c("nodecov.sum.MeanHsig" = "Swell", "absdiff.sum.MeanHsig" = "Swell",
                                                                       "nodecov.sum.TxoccTS" = "Occupancy Tiger Shark", "absdiff.sum.TxoccTS" = "Occupancy Tiger Shark",
                                                                       "nodecov.sum.meanRR.mm." = "Rain", "absdiff.sum.meanRR.mm." = "Rain",
                                                                       "nodecov.sum.MeanDensityTurtle" = "Turtle density", "absdiff.sum.MeanDensityTurtle" = "Turtle density",
                                                                       "nodecov.sum.MeanTotalFreq" = "Human activities", "absdiff.sum.MeanTotalFreq" = "Human activities",
                                                                       "nodecov.sum.rs_moy" = "Turbidity", "absdiff.sum.rs_moy" = "Turbidity",
                                                                       "nodecov.sum.nb_bull" = "Bull Shark abundance", "absdiff.sum.nb_bull" = "Bull Shark abundance",
                                                                       "nodecov.sum.TP_Meme" = "Same sex CRT", "absdiff.sum.TP_Meme" = "Same sex CRT",
                                                                       "nodecov.sum.TP_Diff" = "Opposite sex CRT", "absdiff.sum.TP_Diff" = "Opposite sex CRT",
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

# remove binary networks 8 and 24: remove Jade 16 and Camille 15
ergm.count.tab2<-ergm.count.tab2[-c(1:20,396:406 ),]

df<-ddply(ergm.count.tab2, .(Variables,Term), summarise, mean=mean(Coefficient), se = sd(Coefficient)/sqrt(length(Coefficient)),lowerCI = mean(Coefficient)-qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), upperCI = mean(Coefficient)+qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), n = length(Coefficient))


pd <- position_dodge(0.5) # move them .05 to the left and right

#ggplot(df, aes(x=factor(Variables, level = c('Swell', 'Turbidity', 'Rain','Human activities','Turtle density','Same sex CRT','Opposite sex CRT','Bull Shark abundance', 'Occupancy Bull Shark', 'Tiger Shark abundance','Occupancy Tiger Shark','Geodistance','sum', 'nonzero')), y=mean, color=Term)) + 
ggplot(df, aes(x=factor(Variables, level = c('Geodistance','sum','nonzero','Occupancy Tiger Shark','Tiger Shark abundance','Occupancy Bull Shark','Bull Shark abundance','Opposite sex CRT','Same sex CRT','Human activities','Turtle density','Rain','Turbidity','Swell')), y=mean, color=Term)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0, size=1, position=pd) +
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0, size=0.5, position=pd) +
  geom_point(aes(size=n), alpha=0.5, position=pd)+
  scale_size(breaks=c(1,10,20,30,40), range = c(1,4)) +
  geom_hline(yintercept = 0,linetype = "dotted")+
  scale_color_brewer(palette="Dark2")+
  coord_flip()+
  xlab('Variable')+
  theme_bw()

ergm.count.tab2[which(ergm.count.tab2$Variables=="Opposite sex CRT"),1:5]

range(ergm.count.tab2[which(ergm.count.tab2$Variables=="Opposite sex CRT"),5])



#### Graph for Sex

df_sex<-ddply(ergm.count.tab2, .(Variables,Term,Sex), summarise, mean=mean(Coefficient), se = sd(Coefficient)/sqrt(length(Coefficient)),lowerCI = mean(Coefficient)-qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), upperCI = mean(Coefficient)+qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), n = length(Coefficient))


pd <- position_dodge(0.5) # move them .05 to the left and right
ggplot(df_sex, aes(x=factor(Variables, level = c('Geodistance','sum','nonzero','Occupancy Tiger Shark','Tiger Shark abundance','Occupancy Bull Shark','Bull Shark abundance','Opposite sex CRT','Same sex CRT','Human activities','Turtle density','Rain','Turbidity','Swell')), y=mean, color=Term)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0, size=1, position=pd) +
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0, size=0.5, position=pd) +
  geom_point(aes(size=n), alpha=0.5, position=pd)+
  scale_size(breaks=c(1,10,20,30,40), range = c(1,4)) +
  geom_hline(yintercept = 0,linetype = "dotted")+
  scale_color_brewer(palette="Dark2")+
  coord_flip()+
  facet_grid(~Sex)+
  xlab('Variable')+
  theme_bw()


# Graph for Season

df_season<-ddply(ergm.count.tab2, .(Variables, Term, Season), summarise, mean=mean(Coefficient), se = sd(Coefficient)/sqrt(length(Coefficient)),lowerCI = mean(Coefficient)-qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), upperCI = mean(Coefficient)+qnorm(0.95)*(sd(Coefficient)/sqrt(length(Coefficient))), n = length(Coefficient))

pd <- position_dodge(0.5) # move them .05 to the left and right
ggplot(df_season, aes(x=factor(Variables, level = c('Geodistance','sum','nonzero','Occupancy Tiger Shark','Tiger Shark abundance','Occupancy Bull Shark','Bull Shark abundance','Opposite sex CRT','Same sex CRT','Human activities','Turtle density','Rain','Turbidity','Swell')), y=mean, color=Term)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0, size=1, position=pd) +
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0, size=0.5, position=pd) +
  geom_point(aes(size=n), alpha=0.5, position=pd)+
  scale_size(breaks=c(1,10,20,30,40), range = c(1,4)) +
  geom_hline(yintercept = 0,linetype = "dotted")+
  scale_color_brewer(palette="Dark2")+
  coord_flip()+
  facet_grid(~Season)+
  xlab('Variable')+
  theme_bw()



