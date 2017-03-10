#Simulated dataser of 250 counties if you use your own data, just ignore these lines until library
setwd("c:/CountyBinomial")

#Specified seed for replicability
set.seed(500)

#Number of counties to simulate
nCounties = 250

#Simulates number of species (inflating small numbers, but omiting 0)
N = round(abs(rnorm(nCounties,0,10)))
N[N == 0] = sample(N[N != 0], size = sum(N == 0), replace = T)

#Random P per county from 0 to 1
p = runif(nCounties)

# Generates observations y for N trials and p for each county
y = vector(length=nCounties)
for (i in 1:nCounties)
y[i]=rbinom(size=N[i], prob=p[i], n =1)

#Parameters relating p to X (and its error)
a = -5
b = 10
sigma = 1

#Error (residuals) simulated for each county 
error = rnorm(nCounties,0,sigma)

#X calculated from the logit of p, the error and a and b
logitP= log(p/(1-p))
X = (logitP-error-a)/b
Y = rnorm(nCounties, 10,100)

####END OF SIMULATED DATASET#####



#Load Stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Load your data, N, y and X in this case are vectors without NAs of the same length (number of counties).
#N is the number of species per county (Integer)
#y is the number of endangered species per county (Integer)
#X is the predicting variable (any value) in a matrix format (so 1 or more variables can be used. If only one variable is being used it needs to be a matrix of 1 column)

X = as.matrix(X)
X = cbind(X, Y)

stanData = list( nCounties = length(N),
                 N = N,
                 y = y
                 
)

y_prop = y/N

#This line loads the compiled C++ Stan Model, if all files are in the same working directory, you shouldnt need to compile again, unless you change the .stan file
StanModel = stan_model(file = "OnlyP_countyBinomial_Endemics.stan" )

#Runs the MCMC model                       
FitModel = sampling(StanModel,
data = stanData,              # named list of data
chains = 2,                   # number of Markov chains
warmup = 1000,               # number of warmup iterations per chain
iter = 2000,                 # total number of iterations per chain
cores = 2,                    # number of cores
refresh = 50                 # show progress every 'refresh' iterations
)

#Opens model analytics in Shiny


library(shinystan)
launch_shinystan(FitModel)

