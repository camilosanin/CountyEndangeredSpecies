data{
  int<lower=1> nCounties;      //Number of counties
  int<lower=1> nSpecies;   //Number of species
  int<lower=1> nEnd;   //Number  of endangeredspecies
  int<lower=1> nNotEnd;   //Number of not endangered species
  matrix [nSpecies,nCounties] spOccMat; //Matrix of species occurence
  int<lower=0> endSpp[nEnd];   //which species are enangered
  int<lower=0> notEndSpp[nNotEnd];   //which species are NOT enangered

  
}

transformed data{
  matrix [nEnd,nCounties] endSpOccMat; //Matrix of species occurence
  matrix [nNotEnd,nCounties] notEndSpOccMat; //Matrix of species occurence
  
  endSpOccMat = spOccMat[endSpp,] ; 
  notEndSpOccMat = spOccMat[notEndSpp,] ;

}

parameters{

  vector <upper=0> [nCounties] logP;  

}

transformed parameters {


  vector [nCounties] log1mP;




  log1mP = log1m_exp(logP);
  
}

model
  {
    
  
    //Not Endangered Species (sum of the probabilities that it is not endangered in ANY of the counties)
    
    target +=  notEndSpOccMat * log1mP;
    
    // Endangered Species (see algebra for getting here)
    
    target +=  log1m_exp(endSpOccMat * log1mP);
    

   
    
  }

generated quantities 
  {
  
vector<lower=0, upper=1>[nCounties] p;

p = exp(logP);




}
