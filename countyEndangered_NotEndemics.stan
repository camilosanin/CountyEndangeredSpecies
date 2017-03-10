data{
  int<lower=1> nCounties;      //Number of counties
  int<lower=1> nSpecies;   //Number of species
  int<lower=1> nEnd;   //Number  of endangeredspecies
  int<lower=1> nNotEnd;   //Number of not endangered species
  matrix [nSpecies,nCounties] spOccMat; //Matrix of species occurence
  int<lower=0> endSpp[nEnd];   //which species are enangered
  int<lower=0> notEndSpp[nNotEnd];   //which species are NOT enangered
  int<lower=1> nHyperP; //Number of hyperParameters for p
  int<lower=1> HyperPAssign[nCounties]; //Assignment of counties to hyperPs
  int<lower=1> K;
  matrix [nCounties,K] x_pred; //Matrix of predictors
  
}

transformed data{
  matrix [nEnd,nCounties] endSpOccMat; //Matrix of species occurence
  matrix [nNotEnd,nCounties] notEndSpOccMat; //Matrix of species occurence
  
  endSpOccMat = spOccMat[endSpp,] ; 
  notEndSpOccMat = spOccMat[notEndSpp,] ;

}

parameters{
  real  a;
  vector[K]  b;
  real <lower = 0> sigma;
  vector [nCounties] error;
  vector  [nHyperP] a_cat;
  real <lower = 0> a_sigma;

  

}

transformed parameters {

  vector [nCounties] odds_raw;
  vector<upper = 0> [nCounties] logP;  
  vector<upper = 0> [nCounties] log1mP;



  odds_raw = (x_pred*b)+error+a+a_cat[HyperPAssign];
  logP =  log_inv_logit(odds_raw);
  log1mP = log1m_exp(logP);

}

model
  {
    
    target += cauchy_lpdf(sigma | 0, 1);
    target += normal_lpdf(error | 0, sigma);
    
    //Hyper parameter for P
    
    target += normal_lpdf(a_cat | 0, a_sigma);
    target += cauchy_lpdf(a_sigma | 0, 1);

    //Not Endangered Species (sum of the probabilities that it is not endangered in ANY of the counties)
    
    target +=  notEndSpOccMat * log1mP;
    
    // Endangered Species (see algebra for getting here)
    
    target +=  log1m_exp(endSpOccMat * log1mP);
    

   
    
  }

generated quantities 
  {
  
vector<lower=0, upper=1>[nCounties] sim_p;
vector<lower=0, upper=1>[nCounties] calc_p;
vector<lower=0, upper=1>[nCounties] calc_p_justX;
real  r_sq;
real  r_sq_justX;
vector<lower=0, upper=1>[nCounties] p;


p = exp(logP);




calc_p = inv_logit((x_pred*b)+a+a_cat[HyperPAssign]);
calc_p_justX = inv_logit((x_pred*b)+a) ;

for (cty in 1:nCounties){
  
        sim_p[cty]=inv_logit((x_pred[cty]*b)+normal_rng(0,sigma)+a+a_cat[HyperPAssign[cty]]);

  }

r_sq = 1 - variance(p-calc_p)/variance(p);
r_sq_justX = 1 - variance(p-calc_p_justX)/variance(p);


}
