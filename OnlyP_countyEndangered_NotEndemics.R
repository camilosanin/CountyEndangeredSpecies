#Simulated dataser of 250 counties if you use your own data, just ignore these lines until library
#setwd("c:/CountyBinomial")

#Specified seed for replicability
set.seed(500)

#Number of counties to simulate
nCounties = 500
nSpecies = 2000

#Random P per county from 0 to 1 with a mean of 0.1
#p = rbeta(nCounties, 5,5)
p = runif(nCounties)

#Simulates number of counties a species occurs in (inflating small numbers, but omiting 0)
spCount = round(abs(rnorm(nSpecies,0,3)))
#spCount = rep(1,nSpecies)

spCount[spCount == 0 | spCount > nCounties] = sample(spCount[spCount != 0 & spCount <= nCounties], size = sum(spCount == 0 | spCount > nCounties), replace = T)

#Simulates occurrence matrix (Species by Counties matrix with 1 where they occur 0 where they dont)

occMatrix = matrix(0, ncol = nCounties, nrow = nSpecies)

for (sp in 1:nSpecies){
  
  whichOccurs = sample(1:nCounties, size = spCount[sp])
  
  occMatrix[sp,whichOccurs] = 1 
  
}

#Simulates what species are endangered in at least 1 county

sp = 1
EndSpec = vector(length = nSpecies)

for (sp in 1:nSpecies){
  
  whichOccurs = which(occMatrix[sp,] == 1)
  EndSpec[sp] = max(rbinom(size = 1, n = length(whichOccurs), prob = p[whichOccurs] ))
  
   
  
}





#Parameters relating p to X (and its error)
a = -5
b = 10
sigma = 1

#Error (residuals) simulated for each county 
error = rnorm(nCounties,0,sigma)

#X calculated from the logit of p, the error and a and b
logitP= log(p/(1-p))
X = (logitP-error-a)/b

#Y is a random not correlated variable 
Y = rnorm(nCounties, 10,100)

####END OF SIMULATED DATASET#####



#Load Stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Load your data, N, y and X in this case are vectors without NAs of the same length (number of counties).
#X is the predicting variable (any value) in a matrix format (so 1 or more variables can be used. If only one variable is being used it needs to be a matrix of 1 column)

#Make sure there are objects nCounties and nSpecies with the number of species and counties respectively

#Create x predictors Matrix
X = as.matrix(X)
X = cbind(X, Y)

#Identitity and Unit matrices (I'm sure there is a better way of doing this in Stan (or R) but this should work)

whichEnd = which(EndSpec == 1)
whichNotEnd  = which(EndSpec == 0)

nEnd = length(whichEnd)
nNotEnd = length(whichNotEnd)



stanData = list( nCounties = nCounties,
                 nSpecies = nEnd+nNotEnd,
                 nEnd = nEnd,
                 nNotEnd = nNotEnd,
                 spOccMat = occMatrix,#Species by county Matrix (1 and 0)
                 endSpp = whichEnd , #Which species are endangered (numbers matching spOccMat)
                 notEndSpp = whichNotEnd #Which species are NOT endangered (numbers matching spOccMat)
                 
)


#This line loads the compiled C++ Stan Model, if all files are in the same working directory, you shouldnt need to compile again, unless you change the .stan file
StanModel = stan_model(file = "OnlyP_countyEndangered_NotEndemics.stan" )

#Runs the MCMC model                       
FitModel = sampling(StanModel,
data = stanData,              # named list of data
chains = 2,                   # number of Markov chains
warmup = 1000,               # number of warmup iterations per chain
iter = 2000,                 # total number of iterations per chain
cores = 2,                    # number of cores
refresh = 50                 # show progress every 'refresh' iterations
)

save(FitModel, file= paste("./",format(Sys.time(), "%b%d%Y"),"FitModel_NoEndemics.RData",sep =""))

#Opens model analytics in Shiny


library(shinystan)
launch_shinystan(FitModel)

