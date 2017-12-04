functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  * FORM https://github.com/mbjoseph/CARstan/blob/master/stan/car_sparse.stan
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}


data{
  int<lower=1> nCounties; //Number of cells that have been sampled
  int<lower=1> sampledId [nCounties]; //Id of sampled cells in complete raster
  int<lower=0> nNotSampled; //Number of cells that have not been sampled
  int<lower=1> notSampledId[nNotSampled]; // Id of not sampled cells in complete raster
  int<lower=1> n;      //Number of no NA cells
  int<lower=1> W_n; // Number of adjacent pairs
  int<lower=1> W_sparse[W_n, 2];   // adjacency pairs
  vector<lower=1>[n] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[n] lambda;       // eigenvalues of invsqrtD * W * invsqrtD


  int<lower=1> nSpecies;   //Number of species
  int<lower=1> nEnd;   //Number  of endangeredspecies
  int<lower=1> nNotEnd;   //Number of not endangered species
  matrix [nSpecies,nCounties] spOccMat; //Matrix of species occurence
  int<lower=0> endSpp[nEnd];   //which species are enangered
  int<lower=0> notEndSpp[nNotEnd];   //which species are NOT enangered
  int<lower=1> nHyperP; //Number of hyperParameters for p
  int<lower=1> HyperPAssign[n]; //Assignment of counties to hyperPs
  int<lower=1> K;
  matrix [n,K] x_pred; //Matrix of predictors
  real <lower=0, upper =1> proportionEndSp;
  

  
}

transformed data{
  matrix [nEnd,nCounties] endSpOccMat; //Matrix of species occurence
  matrix [nNotEnd,nCounties] notEndSpOccMat; //Matrix of species occurence
  int <lower = 0, upper =1> boolEndangered[nSpecies];
  real <lower = 0> logloss_naive;

  boolEndangered[notEndSpp] = rep_array(0, nNotEnd);
  boolEndangered[endSpp] = rep_array(1, nEnd);
  endSpOccMat = spOccMat[endSpp,] ; 
  notEndSpOccMat = spOccMat[notEndSpp,] ;
  
  logloss_naive =  -1 * ((proportionEndSp * log(proportionEndSp)) + ((1 - proportionEndSp) * log(1 - proportionEndSp)));


}

parameters{
  real  a;
  vector[K]  b;
  real <lower = 0> sigma;
  vector [nCounties] error;
  vector [n] geo_effect;
  vector  [nHyperP] a_cat;
  real <lower = 0> a_sigma;
  real <lower=0> tau;
  real <lower=0, upper=1> alpha;


  

}

transformed parameters {

  vector [nCounties] odds_raw;
  vector<upper = 0> [nCounties] logP;  
  vector<upper = 0> [nCounties] log1mP;




  odds_raw = (x_pred[sampledId,]*b)+error+a+a_cat[HyperPAssign[sampledId]]+geo_effect[sampledId];
  logP =  log_inv_logit(odds_raw);
  log1mP = log1m_exp(logP);
  
   
}

model
  {
    target += sparse_car_lpdf(geo_effect | tau, alpha, W_sparse, D_sparse, lambda, n, W_n);
    target += normal_lpdf(geo_effect | 0, 1);


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
  

vector<lower=0, upper=1>[n] sim_p;
vector<lower=0, upper=1>[n] calc_p;
vector<lower=0, upper=1>[n] calc_p_notGeo;
vector<lower=0, upper=1>[n] calc_p_justX;
vector<upper=0>[n] calc_p_log;
vector<upper=0>[n] calc_p_notGeo_log;
vector<upper=0>[n] calc_p_justX_log;
real  obs_lh;
real  random_lh;
real<lower = 0>   logloss_obs;
real<lower = 0>   logloss_random;
real<lower = 0>   logloss_calc_p;
real<lower = 0>   logloss_notGeo;
real<lower = 0>   logloss_justX;
real<lower = 0>   logloss_obs_stand;
real<lower = 0>   logloss_random_stand;
real<lower = 0>   logloss_calc_p_stand;
real<lower = 0>   logloss_notGeo_stand;
real<lower = 0>   logloss_justX_stand;
real  r_sq;
real  r_sq_notGeo;
real  r_sq_justX;
real  r_sq_log;
real  r_sq_notGeo_log;
real  r_sq_justX_log;
vector<lower=0, upper=1>[nCounties] p;
vector <upper = 0> [nNotEnd]probSpNotEnd;
vector <upper = 0> [nEnd]probSpEnd;
vector <upper = 0> [nSpecies]allProbSpNotEnd;
vector <upper = 0> [nSpecies]allProbSpEnd;
vector <upper = 0> [nSpecies]allProbSpNotEnd_calc_p;
vector <upper = 0> [nSpecies]allProbSpEnd_calc_p;
vector <upper = 0> [nSpecies]allProbSpNotEnd_notGeo;
vector <upper = 0> [nSpecies]allProbSpEnd_notGeo;
vector <upper = 0> [nSpecies]allProbSpNotEnd_justX;
vector <upper = 0> [nSpecies]allProbSpEnd_justX;
vector <upper = 0> [nSpecies]logloss_obs_sp;
vector <upper = 0> [nSpecies]logloss_random_sp;
vector <upper = 0> [nSpecies]logloss_calc_p_sp;
vector <upper = 0> [nSpecies]logloss_notGeo_sp;
vector <upper = 0> [nSpecies]logloss_justX_sp;

int<lower=1> rndEndSpp[nEnd];   //which species are enangered
int<lower=1> rndNotEndSpp[nNotEnd];   //which species are NOT enangered
int<lower=0, upper = nEnd> i_end;
int<lower=0, upper = nNotEnd> i_notEnd;



probSpNotEnd = notEndSpOccMat * log1mP;
probSpEnd = log1m_exp(endSpOccMat * log1mP);


i_end = 0;
i_notEnd = 0;

allProbSpNotEnd = spOccMat * log1mP;
allProbSpEnd = log1m_exp(spOccMat * log1mP);


for (sp in 1:(nSpecies)){
  if((i_end < nEnd) && bernoulli_rng(proportionEndSp)){

    i_end = i_end+1;
    rndEndSpp [i_end] = sp;

  }else{

    if(i_notEnd < nNotEnd){
    i_notEnd = i_notEnd+1;
    rndNotEndSpp [i_notEnd] = sp;

    }else{
    i_end = i_end+1;
    rndEndSpp [i_end] = sp;

    }

  }


}

logloss_obs_sp[endSpp] = allProbSpEnd[endSpp];
logloss_random_sp[rndEndSpp]= allProbSpEnd[rndEndSpp];

logloss_obs_sp[notEndSpp] = allProbSpNotEnd[notEndSpp];
logloss_random_sp[rndNotEndSpp]= allProbSpNotEnd[rndNotEndSpp];


logloss_obs = (sum(logloss_obs_sp)/nSpecies)*-1;
logloss_random = (sum(logloss_random_sp)/nSpecies)*-1;

random_lh = sum(allProbSpNotEnd[rndNotEndSpp])+sum(allProbSpEnd[rndEndSpp]);
obs_lh =  sum(probSpNotEnd)+sum(probSpEnd);


p = exp(logP);

calc_p_log = log_inv_logit((x_pred*b)+a+a_cat[HyperPAssign]+geo_effect);
calc_p_notGeo_log = log_inv_logit((x_pred*b)+a+a_cat[HyperPAssign]);
calc_p_justX_log = log_inv_logit((x_pred*b)+a) ;

calc_p = exp(calc_p_log);
calc_p_notGeo = exp(calc_p_notGeo_log);
calc_p_justX =  exp(calc_p_justX_log);


for (cty in 1:n){
  
        sim_p[cty]=inv_logit((x_pred[cty]*b)+normal_rng(0,sigma)+a+a_cat[HyperPAssign[cty]]+geo_effect[cty]);

  }
  
  
  
allProbSpNotEnd_calc_p = spOccMat * log1m_exp(calc_p_log[sampledId]);
allProbSpEnd_calc_p = log1m_exp(spOccMat * log1m_exp(calc_p_log[sampledId]));

allProbSpNotEnd_notGeo = spOccMat * log1m_exp(calc_p_notGeo_log[sampledId]);
allProbSpEnd_notGeo = log1m_exp(spOccMat * log1m_exp(calc_p_notGeo_log[sampledId]));

allProbSpNotEnd_justX = spOccMat * log1m_exp(calc_p_justX_log[sampledId]);
allProbSpEnd_justX = log1m_exp(spOccMat * log1m_exp(calc_p_justX_log[sampledId]));

logloss_calc_p_sp[endSpp] = allProbSpEnd_calc_p[endSpp];
logloss_calc_p_sp[notEndSpp] = allProbSpNotEnd_calc_p[notEndSpp];
logloss_notGeo_sp[endSpp] = allProbSpEnd_notGeo[endSpp];
logloss_notGeo_sp[notEndSpp] = allProbSpNotEnd_notGeo[notEndSpp];
logloss_justX_sp[endSpp] = allProbSpEnd_justX[endSpp];
logloss_justX_sp[notEndSpp] = allProbSpNotEnd_justX[notEndSpp];


logloss_calc_p = (sum(logloss_calc_p_sp)/nSpecies)*-1;
logloss_notGeo = (sum(logloss_notGeo_sp)/nSpecies)*-1;
logloss_justX = (sum(logloss_justX_sp)/nSpecies)*-1;



r_sq = 1 - variance(p-calc_p[sampledId])/variance(p);
r_sq_notGeo = 1 - variance(p-calc_p_notGeo[sampledId])/variance(p);
r_sq_justX = 1 - variance(p-calc_p_justX[sampledId])/variance(p);

r_sq_log = 1 - variance(logP-calc_p_log[sampledId])/variance(logP);
r_sq_notGeo_log = 1 - variance(logP-calc_p_notGeo_log[sampledId])/variance(logP);
r_sq_justX_log = 1 - variance(logP-calc_p_justX_log[sampledId])/variance(logP);

logloss_obs_stand = logloss_obs/logloss_naive;
logloss_random_stand = logloss_random/logloss_naive;
logloss_calc_p_stand = logloss_calc_p/logloss_naive;
logloss_notGeo_stand = logloss_notGeo/logloss_naive;
logloss_justX_stand = logloss_justX/logloss_naive;

}
