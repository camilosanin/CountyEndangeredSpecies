#Simulated dataser of 250 counties if you use your own data, just ignore these lines until library
#setwd("c:/CountyBinomial")

source("shapeLopodData.R")
library(rgeos)
library(raster)
library(ggplot2)
library(slam)

#Specified seed for replicability
set.seed(500)

countyShapeFile = shapefile("./County_Shape/cb_2016_us_county_20m.shp")

#Create a feature with P which is autocorrelated in two ways: just resuduals and for X

nCounties = dim(countyShapeFile@data)[1]


#Autocorrelation strenght 
DenseMatrix = gTouches(countyShapeFile,  byid = T, returnDense=T)

aut_coef = 0.90
countyShapeFile@data[,"P_Geo"] = NA

while (sum(is.na(countyShapeFile@data[,"P_Geo"]))!=0){

  focalFeature = sample(rownames(countyShapeFile@data[which(is.na(countyShapeFile@data[,"P_Geo"])),]), size = 1)
  countyShapeFile@data[focalFeature,"P_Geo"] = rbeta(1,2,2)
  AdjFeats = names(DenseMatrix[focalFeature,DenseMatrix[focalFeature,]])
  naAdjFeat = AdjFeats[which(is.na(countyShapeFile@data[AdjFeats,"P_Geo"]))]

  while(length(naAdjFeat)>0){

    preP = countyShapeFile@data[focalFeature,"P_Geo"]
    countyShapeFile@data[naAdjFeat,"P_Geo"] = preP*aut_coef + runif(length(naAdjFeat))*(1-aut_coef)
    focalFeature = sample(naAdjFeat, size = 1)
    AdjFeats = names(DenseMatrix[focalFeature,DenseMatrix[focalFeature,]])
    naAdjFeat = AdjFeats[which(is.na(countyShapeFile@data[AdjFeats,"P_Geo"]))]

    print(sum(is.na(countyShapeFile@data[,"P_Geo"])))

  }
  message("JUMP!")

}

seed(500)
aut_coef = 0.95
countyShapeFile@data[,"P_X"] = NA

while (sum(is.na(countyShapeFile@data[,"P_X"]))!=0){
  
  focalFeature = sample(rownames(countyShapeFile@data[which(is.na(countyShapeFile@data[,"P_X"])),]), size = 1)
  countyShapeFile@data[focalFeature,"P_X"] = rbeta(1,0.5,0.5)
  AdjFeats = names(DenseMatrix[focalFeature,DenseMatrix[focalFeature,]])
  naAdjFeat = AdjFeats[which(is.na(countyShapeFile@data[AdjFeats,"P_X"]))]
  
  while(length(naAdjFeat)>0){
    
    preP = countyShapeFile@data[focalFeature,"P_X"]
    countyShapeFile@data[naAdjFeat,"P_X"] = preP*aut_coef + runif(length(naAdjFeat))*(1-aut_coef)
    focalFeature = sample(naAdjFeat, size = 1)
    AdjFeats = names(DenseMatrix[focalFeature,DenseMatrix[focalFeature,]])
    naAdjFeat = AdjFeats[which(is.na(countyShapeFile@data[AdjFeats,"P_X"]))]
    
    print(sum(is.na(countyShapeFile@data[,"P_X"])))
    
  }
  message("JUMP!")
  
}
# Contribution from residuals
countyShapeFile@data[,"P_Res"] = runif(nCounties)

#Contrubution from State
countyShapeFile@data[,"P_State"] =  NA

for(st in unique(countyShapeFile@data$STATEFP) ){
  
  countyShapeFile@data[countyShapeFile@data$STATEFP==st,"P_State"] = runif(1)
  
  
}

#All effects combined in P

countyShapeFile@data[,"P"] = countyShapeFile@data[,"P_X"]*0.55+countyShapeFile@data[,"P_Geo"]*0.20+countyShapeFile@data[,"P_State"]*0.15+countyShapeFile@data[,"P_Res"]*0.1

#spplot(countyShapeFile,  col = "transparent", zcol = c("P_Geo", "P_X", "P_State", "P"))
spplot(countyShapeFile,  col = "transparent", zcol = c("P"))

#Number of species to simulate
nSpecies = 2500


#Simulates number of counties a species occurs in (inflating small numbers, but omiting 0)
spCount = round(abs(rnorm(nSpecies,0,5)))
#spCount = rep(1,nSpecies)

spCount[spCount == 0 | spCount > nCounties] = sample(spCount[spCount != 0 & spCount <= nCounties], size = sum(spCount == 0 | spCount > nCounties), replace = T)

#Simulates occurrence matrix (Species by Counties matrix with 1 where they occur 0 where they dont)

sp = 1
SppCols = vector()

for (sp in 1:nSpecies){
  colNAme = paste("Sp", sp, sep="_")
  SppCols = c(SppCols,colNAme)
  countyShapeFile@data[,colNAme] = 0
  nCountSp = spCount[sp]
  
  focalCounty = sample(rownames(countyShapeFile@data), size = 1)
  countyShapeFile@data[focalCounty,colNAme] = 1
  nCountSp = nCountSp - 1
  
  
  while(nCountSp > 0){
    AdjFeats = names(DenseMatrix[focalCounty,DenseMatrix[focalCounty,]])
    naAdjFeat = AdjFeats[which(countyShapeFile@data[AdjFeats,colNAme]==0)]
    
    if (length(naAdjFeat)>0){
    
    if(length(naAdjFeat)>nCountSp){
      naAdjFeat = sample(naAdjFeat,nCountSp)
    }
    
    countyShapeFile@data[naAdjFeat,colNAme] = 1
    nCountSp = nCountSp - length(naAdjFeat)
    focalCounty = sample(naAdjFeat,1)
    
    } else {
      
      focalCounty = sample(rownames(countyShapeFile@data[countyShapeFile@data[,colNAme]==0,]), size = 1)
      countyShapeFile@data[focalCounty,colNAme] = 1
      nCountSp = nCountSp - 1
      message("JUMP!")
      
    }
    
    }
  
print(paste(colNAme, "of",nSpecies ))
}

countyShapeFile@data[,"Richness"] = rowSums(countyShapeFile@data[,SppCols])
spplot(countyShapeFile, zcol="Richness", col="transparent")


which(countyShapeFile@data[,"Richness"] ==0)

#Simulates what species are endangered in at least 1 county

sp = 1
EndSpec = vector(length = nSpecies)
names(EndSpec) = SppCols


for (sp in SppCols){
  
  whichOccurs = which(countyShapeFile@data[,sp] == 1)
  EndSpec[sp] = max(rbinom(size = 1, n = length(whichOccurs), prob = countyShapeFile@data[whichOccurs,"P"] ))
  
  
}

hist(spCount)
hist(EndSpec)


#Parameters relating p to X (and its error)
a = -2.5
b = 5



#X calculated from the logit of p, the error and a and b
logitP= log(countyShapeFile@data[,"P_X"]/(1-countyShapeFile@data[,"P_X"]))
countyShapeFile@data[,"X"] = (logitP-a)/b

#Y is a random not correlated variable 
countyShapeFile@data[,"Y"] = rnorm(nCounties, 10,100)


#spplot(countyShapeFile, col ="transparent", zcol = "X")
#spplot(countyShapeFile, col ="transparent", zcol = "Y")


####END OF SIMULATED DATASET#####


AdMatrixData = shapeLopodData(countyShapeFile, fieldN = "Richness", fieldY = "Richness", keepFields = TRUE, Adjacency = TRUE)


#Load Stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Load your data, N, y and X in this case are vectors without NAs of the same length (number of counties).
#X is the predicting variable (any value) in a matrix format (so 1 or more variables can be used. If only one variable is being used it needs to be a matrix of 1 column)

#Make sure there are objects nCounties and nSpecies with the number of species and counties respectively

AllIncludedFeatures = rbind(AdMatrixData$geoInfo$sampledId,AdMatrixData$geoInfo$notSampledId)
AllIncludedFeatures = AllIncludedFeatures[order(AllIncludedFeatures$cellStan),]

#Create x predictors Matrix
X = as.matrix(cbind(AdMatrixData$geoDataObject@data[AllIncludedFeatures$featureShape,"X"], AdMatrixData$geoDataObject@data[AllIncludedFeatures$featureShape,"Y"]))


whichEnd = which(EndSpec == 1)
whichNotEnd  = which(EndSpec == 0)

nEnd = length(whichEnd)
nNotEnd = length(whichNotEnd)

#Field in shapefile with state
nHyperP = length(unique(AdMatrixData$geoDataObject@data$STATEFP))
HyperPAssign = as.numeric(factor(as.character(AdMatrixData$geoDataObject@data$STATEFP[AllIncludedFeatures$featureShape])))

spOccMat = t(AdMatrixData$geoDataObject@data[AdMatrixData$geoInfo$sampledId$featureShape,SppCols])

#Names of the list elements should be the same (left side of the =)
stanData = list( nCounties = length(AdMatrixData$geoInfo$sampledId$featureShape),
                 sampledId = AdMatrixData$geoInfo$sampledId$cellStan,
                 nNotSampled = length(AdMatrixData$geoInfo$notSampledId$featureShape),
                 notSampledId = AdMatrixData$geoInfo$notSampledId$cellStan,
                 n = length(AdMatrixData$geoInfo$sampledId$featureShape) + length(AdMatrixData$geoInfo$notSampledId$featureShape) ,
                 W_n = dim(AdMatrixData$geoInfo$W_sparse)[1] ,
                 W_sparse = AdMatrixData$geoInfo$W_sparse ,
                 D_sparse = AdMatrixData$geoInfo$D_sparse ,
                 lambda = AdMatrixData$geoInfo$lambda_sparse  ,
                 nSpecies = nEnd+nNotEnd,
                 nEnd = nEnd,
                 nNotEnd = nNotEnd,
                 spOccMat = spOccMat,#Species by county Matrix (1 and 0)
                 endSpp = whichEnd , #Which species are endangered (numbers matching spOccMat)
                 notEndSpp = whichNotEnd, #Which species are NOT endangered (numbers matching spOccMat)
                 nHyperP = nHyperP, #Number of categories (States)
                 HyperPAssign = HyperPAssign, #Vector of length counties assigning it to a state (numeric)
                 K = dim(X)[2],
                 x_pred = X
                 
)


#This line loads the compiled C++ Stan Model, if all files are in the same working directory, you shouldnt need to compile again, unless you change the .stan file
StanModel = stan_model(file = "countyEndangered_NotEndemics_CAR.stan" )

#Runs the MCMC model                       
FitModel = sampling(StanModel,
data = stanData,              # named list of data
chains = 2,                   # number of Markov chains
warmup = 1000,               # number of warmup iterations per chain
iter = 2000,                 # total number of iterations per chain (includes warm-up)
cores = 2,                    # number of cores
refresh = 10                 # show progress every 'refresh' iterations
)

#save(FitModel, file= paste("./",format(Sys.time(), "%b%d%Y"),"FitModel_NoEndemics.RData",sep =""))

#Opens model analytics in Shiny

resultsShape = AdMatrixData$geoDataObject
  
library(shinystan)
#launch_shinystan(FitModel)

#Calls coeff for predictors in matrix X and intersect
summary(FitModel, pars = c("a",paste("b[",1:dim(X)[2],"]", sep="")))$summary #[,"50%"]


#Calls coeff for states effects in matrix X
summary(FitModel, pars = paste("a_cat[",1:nHyperP,"]", sep=""))$summary[,"50%"]

resultsShape@data[AllIncludedFeatures$featureShape,"state_coeff"] = summary(FitModel, pars =paste("a_cat[",1:nHyperP,"]", sep=""))$summary[paste("a_cat[",HyperPAssign,"]", sep=""),"50%"]

plot(resultsShape@data[,c("P_State","state_coeff")])
abline(h=0,lty="dashed", lwd= 2, col="darkred")

plot(resultsShape@data[,c("P_Geo","state_coeff")])
abline(h=0,lty="dashed", lwd= 2, col="darkred")



#Calls r squareds 
summary(FitModel, pars = c("r_sq","r_sq_justX","r_sq_notGeo" ))$summary #[,"50%"]


#Calls p per county 

resultsShape@data[AdMatrixData$geoInfo$sampledId$featureShape,"p_est"] = summary(FitModel, pars = paste("p[",1:length(AdMatrixData$geoInfo$sampledId$featureShape),"]", sep=""))$summary[,"50%"]
plot(resultsShape@data[AdMatrixData$geoInfo$sampledId$featureShape,c("P","p_est")])
abline(a = 0,b = 1,lty="dashed", lwd= 2, col="darkred")


#Calls p_calc per county (p calculated from X, a and state)
resultsShape@data[AllIncludedFeatures$featureShape,"p_calc"] = summary(FitModel, pars = paste("calc_p[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$summary[,"50%"]
plot(resultsShape@data[AllIncludedFeatures$featureShape,c("P","p_calc")])
abline(a = 0,b = 1,lty="dashed", lwd= 2, col="darkred")


#Calls p_calc per county (p calculated from X, a NOT state)
resultsShape@data[AllIncludedFeatures$featureShape,"p_calc_justX"] = summary(FitModel, pars = paste("calc_p_justX[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$summary[,"50%"]
plot(resultsShape@data[AllIncludedFeatures$featureShape,c("P","p_calc_justX")])
abline(a = 0,b = 1,lty="dashed", lwd= 2, col="darkred")

#Calls p_calc per county (p calculated from X, a and state NOT Geo)
resultsShape@data[AllIncludedFeatures$featureShape,"p_calc_notGeo"] = summary(FitModel, pars = paste("calc_p_notGeo[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$summary[,"50%"]
plot(resultsShape@data[AllIncludedFeatures$featureShape,c("P","p_calc_notGeo")])
abline(a = 0,b = 1,lty="dashed", lwd= 2, col="darkred")

#Calls geo effect
resultsShape@data[AllIncludedFeatures$featureShape,"geo_effect"] = summary(FitModel, pars = paste("geo_effect[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$summary[,"50%"]
plot(resultsShape@data[AllIncludedFeatures$featureShape,c("P_Geo","geo_effect")])
abline(h=0,lty="dashed", lwd= 2, col="darkred")


#The residuals can be calculated by p - calc_p

stan_dens(FitModel,pars = c("r_sq","r_sq_justX"),separate_chains = T )

stan_dens(FitModel,pars = c("p[1]","p[50]"),separate_chains = T )

shapefile(resultsShape, "./County_Shape/countyOutput_CAR.shp")


