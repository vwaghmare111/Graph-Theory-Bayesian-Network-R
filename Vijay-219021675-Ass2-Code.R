library(igraph)
library(ggm )
library(Rgraphviz)
library(RBGL)
library(gRbase)
library(gRain)
library (bnlearn)
#install.packages("Rgraphviz")
#source("http://bioconductor.org/biocLite.R")
BiocManager::install(c("graph", "RBGL", "Rgraphviz"))

#############################################################################################################

#Question 1.5
dag <- DAG(W~H,N~H+W,S~N,R~H+N,M~H)
plotGraph (dag , nodesize =20 , tcltk =FALSE , vc=" white ")

dSep (dag , first ="M", second ="S", cond = NULL)
dSep (dag , first ="W", second ="R", cond =c("N","H"))
dSep (dag , first =c("R","S"), second ="W", cond ="H")

#############################################################################################################

#Question 2
#2.1.a Obtaining belief network
A<-cptable(~season,values =c(0.3,0.7),levels=c('wet','dry'))
B<-cptable(~rate,values=c(0.2,0.8),levels=c('high','low'))
C.AB<-cptable(~fish|season:rate,values=c(0.4,0.6,0.5,0.5,0.6,0.4,0.3,0.7),levels=c('bass','cod'))
D.C<-cptable(~colour|fish, values=c(0.2,0.4,0.4,0.5,0.3,0.2),levels = c('light','medium','dark'))
E.C<-cptable(~size|fish,values=c(0.6,0.4,0.4,0.6),levels=c('wide','thin'))


#Compile list of conditional probability tables and create the network:
prob_list <- compileCPT(list(A,B,C.AB,D.C,E.C))
prob_list

#Plotting belief network
network <- grain(prob_list)
summary(network)
#plot DAG
plot(network$dag)

#Marginal probabilities of all nodes
querygrain(network)

#2.1.b Probability tables
prob_list$season
prob_list$rate
prob_list$fish
prob_list$colour
prob_list$size

#2.2 Probabilities
a <- setEvidence(network, evidence=list(rate="low"))
querygrain(a, nodes=c("size"))

b <- setEvidence(network, evidence=list(colour="dark",season="dry"))
querygrain(b,nodes=c("fish"))

c<-querygrain(network, nodes=c("colour","fish"), type="joint" )
c

d<-querygrain(network, nodes=c("fish"), type="marginal")
d

#############################################################################################################

#Question 4

#4.1
#Downloading hailfinder dataset
original_dataset<-hailfinder

#Viewing the dataset
#View(original_dataset)

#Checing summary of the complete datset
dim(original_dataset)
str(hailfinder)
summary(hailfinder)



#creating and plotting the true network structure
modelstring = paste0("[N07muVerMo][SubjVertMo][QGVertMotion][SatContMoist][RaoContMoist]",
                     "[VISCloudCov][IRCloudCover][AMInstabMt][WndHodograph][MorningBound][LoLevMoistAd][Date]",
                     "[MorningCIN][LIfr12ZDENSd][AMDewptCalPl][LatestCIN][LLIW]",
                     "[CombVerMo|N07muVerMo:SubjVertMo:QGVertMotion][CombMoisture|SatContMoist:RaoContMoist]",
                     "[CombClouds|VISCloudCov:IRCloudCover][Scenario|Date][CurPropConv|LatestCIN:LLIW]",
                     "[AreaMesoALS|CombVerMo][ScenRelAMCIN|Scenario][ScenRelAMIns|Scenario][ScenRel34|Scenario]",
                     "[ScnRelPlFcst|Scenario][Dewpoints|Scenario][LowLLapse|Scenario][MeanRH|Scenario]",
                     "[MidLLapse|Scenario][MvmtFeatures|Scenario][RHRatio|Scenario][SfcWndShfDis|Scenario]",
                     "[SynForcng|Scenario][TempDis|Scenario][WindAloft|Scenario][WindFieldMt|Scenario]",
                     "[WindFieldPln|Scenario][AreaMoDryAir|AreaMesoALS:CombMoisture]",
                     "[AMCINInScen|ScenRelAMCIN:MorningCIN][AMInsWliScen|ScenRelAMIns:LIfr12ZDENSd:AMDewptCalPl]",
                     "[CldShadeOth|AreaMesoALS:AreaMoDryAir:CombClouds][InsInMt|CldShadeOth:AMInstabMt]",
                     "[OutflowFrMt|InsInMt:WndHodograph][CldShadeConv|InsInMt:WndHodograph][MountainFcst|InsInMt]",
                     "[Boundaries|WndHodograph:OutflowFrMt:MorningBound][N34StarFcst|ScenRel34:PlainsFcst]",
                     "[CompPlFcst|AreaMesoALS:CldShadeOth:Boundaries:CldShadeConv][CapChange|CompPlFcst]",
                     "[InsChange|CompPlFcst:LoLevMoistAd][CapInScen|CapChange:AMCINInScen]",
                     "[InsSclInScen|InsChange:AMInsWliScen][R5Fcst|MountainFcst:N34StarFcst]",
                     "[PlainsFcst|CapInScen:InsSclInScen:CurPropConv:ScnRelPlFcst]")
true_network = model2network(modelstring)
par(mfrow = c(1,1))
graphviz.plot(true_network)


#Creating 3 subsets with lenth 100, 1000, and 10,000
first_100<-original_dataset[1:100,]
first_1000<-original_dataset[1:1000,]
first_10000<-original_dataset[1:10000,]

#Checing dimension of each subset
dim(first_100)
dim(first_1000)
dim(first_10000)

#Using hill climbing algorithm and BIC score to obtain and plot the network for each subsets
bicnet100<-hc(first_100,score="bic")
bicnet1000<-hc(first_1000,score="bic")
bicnet10000<-hc(first_10000,score="bic")


#Using hill climbing algorithm and BDE score to obtain and plot the network for each subsets
bdenet100<-hc(first_100,score="bde")
bdenet1000<-hc(first_1000,score="bde")
bdenet10000<-hc(first_10000,score="bde")


#Summary of each fit above
bicnet100
bicnet1000
bicnet10000
bicnet20000

bdenet100
bdenet1000
bdenet10000
bdenet20000

#Score of the fits
score_bicnet100<-score(bicnet100,first_100, type="bic")
score_bicnet100
score_bicnet1000<-score(bicnet1000,first_1000, type="bic")
score_bicnet1000
score_bicnet10000<-score(bicnet10000,first_10000, type="bic")
score_bicnet10000


score_bdenet100<-score(bdenet100,first_100, type="bde")
score_bdenet100
score_bdenet1000<-score(bdenet1000,first_1000, type="bde")
score_bdenet1000
score_bdenet10000<-score(bdenet10000,first_10000, type="bde")
score_bdenet10000


#Plott network of each subset
plot(bicnet100, main="First 100_BIC")
plot(bicnet1000,main="First 1000_BIC")
plot(bicnet10000,main="First 10000_BIC")

plot(bdenet100, main="First 100_BDe")
plot(bdenet1000,main="First 1000_BDe")
plot(bdenet10000,main="First 10000_BDe")

#4.3

#Fitting model on entire dataset
#Using hill climbing algorithm and BIC/BDE score to obtain and plot the network for each subsets
bicnet20000<-hc(original_dataset,score="bic")
bdenet20000<-hc(original_dataset,score="bde")

#Summary of each fit above
bicnet20000
bdenet20000

#Score of the fits
score_bicnet20000<-score(bicnet20000,original_dataset, type="bic")
score_bicnet20000
score_bdenet20000<-score(bdenet20000,original_dataset, type="bde")
score_bdenet20000

#Plott network of each subset
plot(bicnet20000,main="Complete_Data_BICe")
plot(bdenet20000,main="Complete_Data_BDe")

#Comparing complete dataset network with true network

#Checking if two networks has same structure
all.equal(true_network,bicnet20000)
all.equal(true_network,bdenet20000)

#Comparing two networks
unlist(compare(true_network,bicnet20000))
unlist(compare(true_network,bdenet20000))

#Structural distance
hamming(true_network,bicnet20000)
hamming(true_network,bdenet20000)

#Graphical comparison
graphviz.compare(true_network,bicnet20000)
graphviz.compare(true_network,bdenet20000)

#Fitting data to the using (BIC)
set.seed(1)
fittedParams = bn.fit(bicnet20000, original_dataset)
fittedParams$CombClouds

#Finding probability of P(CombClouds ="Cloudy" | MeanRH = "VeryMoist", IRCloudCover ="Cloudy"):
cpquery(fittedParams, event = (CombClouds=="Cloudy"), 
        evidence = ((MeanRH=="VeryMoist") & (IRCloudCover =="Cloudy")))


