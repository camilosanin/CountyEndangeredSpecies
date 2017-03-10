data{
  int<lower=1> nCounties;      //Number of counties
  int<lower=1> N[nCounties];   //Number of species
  int<lower=0> y[nCounties];   //Number of endangered species
  int<lower=1> K;
  matrix [nCounties,K] x_pred; //Matrix of predictors
}

transformed data{


}

parameters{
  real  a;
  vector[K]  b;
  real <lower = 0> sigma;
  vector [nCounties] error;
  

}

transformed parameters {

  vector<lower=0, upper=1> [nCounties] p;


  p =  inv_logit((x_pred*b)+error+a);




}

model
  {
    target += cauchy_lpdf(sigma | 0, 1);
    target += normal_lpdf(error | 0, sigma);

    
    target += binomial_lpmf(y | N, p);


    
    
  }

generated quantities 
  {
int<lower=0> sim_y[nCounties]; //Simulated Sampling
vector<lower=0, upper=1>[nCounties] sim_p;


for (cty in 1:nCounties){
  
        sim_p[cty]=inv_logit((x_pred[cty]*b)+normal_rng(0,sigma)+a);

         sim_y[cty]=binomial_rng(N[cty],sim_p[cty]);
      
      

  }


}
