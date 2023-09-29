functions {
	real[] logistic_model(real y0, real fitness_clades_pre_vacc, real fitness_clades_post_vacc, int nb_years, int vaccine_introduction, int ACV_introduction){
		real res[nb_years];
		real tmp;
    tmp = y0;
	  if(tmp < 1E-30) {tmp = 1E-30;}
	  if(tmp >= 1-1E-30) {tmp = 1-1E-30;}
	  for(i in 1:nb_years){
	      res[i] = tmp;
	      if(i <= vaccine_introduction && i<= ACV_introduction){
	        res[i] = tmp;
	      }
	      if(i > vaccine_introduction && i <=  ACV_introduction){
	        if(fitness_clades_pre_vacc<=0){ // 1 - usual logistic
	          res[i] = 1 - 1/(1+(1/(1-tmp)-1)*exp(fitness_clades_pre_vacc*(1)));
	        }
	        if(fitness_clades_pre_vacc > 0){ // usual logistic
	          res[i] = 1/(1+(1/tmp-1)*exp(-fitness_clades_pre_vacc*(1)));
	        }
  	    }
  	    if(i > vaccine_introduction && i > ACV_introduction){
  	      if(fitness_clades_post_vacc<=0){ // 1 - usual logistic
	          res[i] = 1 - 1/(1+(1/(1-tmp)-1)*exp(fitness_clades_post_vacc*(1)));
	        }
	        if(fitness_clades_post_vacc > 0){ // usual logistic
	          res[i] = 1/(1+(1/tmp-1)*exp(-fitness_clades_post_vacc*(1)));
	        }
  	    }
    	  tmp = res[i];
	      if(res[i] <= 1E-30) {res[i] = 1E-30;}
	      if(res[i] >= 1) {res[i] = 1-1E-30;}
	  }
  	return (res);
	}
}

data {
  // Shape data
	int nb_GPSC; // number of GPSC
  int nb_sero; // number of serotypes
  int dim_data; // dimension data
	int nb_years; // length of the data 
  
  // Data
	int data_genotype_non_ref[dim_data-1, nb_years]; // counts
	int data_genotype_ref[nb_years]; // counts 
	int data_total_number[nb_years]; // reference for each year // real because of the divison in lp for the freq
	int non_zero_country_year[nb_years];
	int non_zero_country_year_genotype[dim_data-1, nb_years];
	int number_zeros_country_year;
	int number_zeros_country_year_genotype;
	
  // ACV implementation
	int yearIntroduction; // dates of introdution for each country
	
	// F0 estimation
	int yearF0;
	
	// Fitness vector
  int R_every_pre_vacc;
  int R_every_post_vacc;
  int number_R_pre_vacc;
  int number_R_post_vacc;
  
  // Fitness values: for now constant
  real<lower=-5, upper=100> fitness_genotypes_pre_vacc[dim_data-1];
	real<lower=-5, upper=100> fitness_genotypes_post_vacc[dim_data-1];
}

parameters {
	simplex[dim_data] f;
	real<lower=1E-30, upper=1-1E-30> alpha;
}

transformed parameters{
  real<lower=-5, upper=100> fitness_genotypes_vector_pre_vacc[dim_data-1];
	real<lower=-5, upper=100> fitness_genotypes_vector_post_vacc[dim_data-1];
	
	simplex[dim_data] f0;
		
	real<lower=0, upper=1> pred_paired_freq[dim_data-1, nb_years] ; // pred frequencies
	real pred_paired_freq_tmp[dim_data, nb_years] ; // pred numbers
	real pred_absolute_freq[dim_data, nb_years] ; // pred numbers
	
	real pred_number_non_ref[dim_data-1, nb_years] ; // pred numbers non ref
	real pred_number_ref[ nb_years] ; // pred numbers ref
  
  f0 = f;

	// construct fitness vector
	fitness_genotypes_vector_pre_vacc = fitness_genotypes_pre_vacc;
	fitness_genotypes_vector_post_vacc = fitness_genotypes_post_vacc;
  
  // compute pred fequencies, by GSPC-ser
  for (j in 1:(dim_data-1)){
    if(fitness_genotypes_vector_pre_vacc[j] < 20){
      pred_paired_freq[j,] = logistic_model(f0[j]/(f0[j]+f0[dim_data]), fitness_genotypes_vector_pre_vacc[j], fitness_genotypes_vector_post_vacc[j], nb_years, yearF0, yearIntroduction);
    }
    if(fitness_genotypes_vector_pre_vacc[j] > 20){ // Reference serotype
      pred_paired_freq[j,] = logistic_model(f0[j]/(f0[j]+f0[dim_data]), 0, 0, nb_years, yearF0, yearIntroduction);
    }
  }
  for(l in 1:nb_years){
	    for(j in 1:(dim_data-1)){
          pred_paired_freq_tmp[j,l] = (pred_paired_freq[j,l])*10000/(1-pred_paired_freq[j,l]);
          if(pred_paired_freq_tmp[j,l] < 0) pred_paired_freq_tmp[j,l] = 0;
  	  }
  	  pred_paired_freq_tmp[dim_data,l] = 10000;
  	  for(j in 1:(dim_data)){
  	    pred_absolute_freq[j,l] = 0;
  	    pred_absolute_freq[j,l] = pred_paired_freq_tmp[j,l]/sum(pred_paired_freq_tmp[,l]);
  	    if(pred_absolute_freq[j,l]<0) pred_absolute_freq[j,l] = 0;
  	  }
  	  
  	  for(j in 1:(dim_data-1)){
  	    pred_number_non_ref[j,l] = pred_paired_freq[j,l]*(data_genotype_non_ref[j,l] + data_genotype_ref[l]); 
  	    if (pred_number_non_ref[j,l] <1E-30) {
  		  	pred_number_non_ref[j,l] = 1E-30;
  		  }
      }
	    // Compute ref clade numbers
  	  pred_number_ref[l] = 0;
  	  pred_number_ref[l] = data_total_number[l] - sum(pred_number_non_ref[,l]);
  	  if (pred_number_ref[l] <1E-30) {
  		  pred_number_ref[l] = 1E-30;
  		}
  }
}

model {
  for(l in 1:yearIntroduction){
	    for (j in 1:(dim_data-1)){
	      if(non_zero_country_year_genotype[j,l]) {
	 		    data_genotype_non_ref[j,l] ~ poisson(pred_number_non_ref[j,l]);
	      }
  	  }
  }
}

generated quantities {
  vector[(dim_data-1)*nb_years-number_zeros_country_year_genotype] log_lik;
  int n;
  n = 1;
  // for(l in 1:yearIntroduction){
  for(l in 1:nb_years){
    for (j in 1:(dim_data-1)){
      if(non_zero_country_year_genotype[j,l]) {
	      log_lik[n] = 1E-30;
 		    log_lik[n] = poisson_lpmf(data_genotype_non_ref[j,l]| pred_number_non_ref[j,l]);
 		    n=n+1;
      }
    }
	}
}

