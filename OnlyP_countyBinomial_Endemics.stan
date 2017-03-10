data{
  int<lower=1> nCounties;      //Number of counties
  int<lower=1> N[nCounties];   //Number of species
  int<lower=0> y[nCounties];   //Number of endangered species

}

transformed data{


}

parameters{

  vector<lower=0, upper=1> [nCounties] p;


}

transformed parameters {


}

model
  {
   
    target += binomial_lpmf(y | N, p);


    
    
  }

generated quantities 
  {
int<lower=0> sim_y[nCounties]; //Simulated Sampling


for (cty in 1:nCounties){
  

         sim_y[cty]=binomial_rng(N[cty],p[cty]);
      
      

  }


}
