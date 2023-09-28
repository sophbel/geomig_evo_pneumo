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
	int nb_genotypes_amr; // number of clades
	int nb_genotypes_vstat; // number of clades
	int nb_genotypes_total; // number of clades
	int nb_years; // length of the data 
	int nb_countries; // number of countries 
  
  // Data AMR
	int data_genotype_non_ref_amr[nb_genotypes_amr-1, nb_years, nb_countries]; // counts
	int data_genotype_ref_amr[nb_years, nb_countries]; // counts 
	int data_total_number_amr[nb_years, nb_countries]; // reference for each year // real because of the divison in lp for the freq
	int non_zero_country_year_amr[nb_years, nb_countries];
	int non_zero_country_year_genotype_amr[nb_genotypes_amr-1, nb_years, nb_countries];
	int number_zeros_country_year_amr;
	int number_zeros_country_year_genotype_amr;
	
	// Data Vaccine Status
	int data_genotype_non_ref_vstat[nb_genotypes_vstat-1, nb_years, nb_countries]; // counts
	int data_genotype_ref_vstat[nb_years, nb_countries]; // counts 
	int data_total_number_vstat[nb_years, nb_countries]; // reference for each year // real because of the divison in lp for the freq
	int non_zero_country_year_vstat[nb_years, nb_countries];
	int non_zero_country_year_genotype_vstat[nb_genotypes_vstat-1, nb_years, nb_countries];
	int number_zeros_country_year_vstat;
	int number_zeros_country_year_genotype_vstat;
	
		// Data Vaccine Status and AMR
	int data_genotype_total[nb_genotypes_total, nb_years, nb_countries]; // counts
	int data_total_number_total[nb_years, nb_countries]; // reference for each year // real because of the divison in lp for the freq
	int non_zero_country_year_total[nb_years, nb_countries];
	int non_zero_country_year_genotype_total[nb_genotypes_total, nb_years, nb_countries];
	int number_zeros_country_year_total;
	int number_zeros_country_year_genotype_total;
	
  // ACV implementation
	int yearIntroduction[nb_countries]; // dates of introdution for each country
	
	// F0 estimation
	int yearF0[nb_countries];

}

parameters {
  
  simplex[nb_genotypes_total] f0_total[nb_countries];
  
  //AMR - Single parameter
	// real<lower=-5, upper=5> fitness_genotypes_pre_vacc_amr[nb_genotypes_amr-1, 1];
	// real<lower=-5, upper=5> fitness_genotypes_post_vacc_amr[nb_genotypes_amr-1, 1];
	
	//AMR - Two parameter
	real<lower=-5, upper=5> fitness_genotypes_pre_vacc_amr[nb_genotypes_amr, 1];
	real<lower=-5, upper=5> fitness_genotypes_post_vacc_amr[nb_genotypes_amr, 1];
	//
	// simplex[nb_genotypes_amr] f0_amr[nb_countries];

  //Vaccine Status
  real<lower=-5, upper=5> fitness_genotypes_pre_vacc_vstat[nb_genotypes_vstat-1, 1];
	real<lower=-5, upper=5> fitness_genotypes_post_vacc_vstat[nb_genotypes_vstat-1, 1];
	
	// simplex[nb_genotypes_vstat] f0_vstat[nb_countries];
	
  // real<lower=0.25, upper=10> overdispersion_inv_total;
}

transformed parameters{
  
  ////////////////////////AMR PARAMETERS
  ///single parameter estimate for amr
  // real<lower=-5, upper=5> fitness_genotypes_vector_pre_vacc_amr[nb_genotypes_amr-1, 1];
	// real<lower=-5, upper=5> fitness_genotypes_vector_post_vacc_amr[nb_genotypes_amr-1, 1];
	
	///allowing two parameter estimates for amr
	real<lower=-5, upper=5> fitness_genotypes_vector_pre_vacc_amr[nb_genotypes_amr, 1];
	real<lower=-5, upper=5> fitness_genotypes_vector_post_vacc_amr[nb_genotypes_amr, 1];
	
	
	///single parameter estimate for amr
	// real<lower=0, upper=1> pred_paired_freq_amr[nb_countries,nb_genotypes_amr-1, nb_years] ; // pred frequencies
	// 
	// real pred_paired_freq_tmp_amr[nb_countries,nb_genotypes_amr, nb_years] ; // pred numbers
	// real pred_absolute_freq_amr[nb_countries,nb_genotypes_amr, nb_years] ; // pred numbers
	// 
	// real pred_number_non_ref_amr[nb_countries,nb_genotypes_amr-1, nb_years] ; // pred numbers non ref
	// real pred_number_ref_amr[nb_countries, nb_years] ; // pred numbers ref
		// simplex[nb_genotypes_amr] f0_amr[nb_countries];

	
 ///adding extra parameter for vaccine status resistance
	real<lower=0, upper=1> pred_paired_freq_amr[nb_countries,nb_genotypes_amr-1, nb_years,nb_genotypes_vstat] ; // pred frequencies
	
	real pred_paired_freq_tmp_amr[nb_countries,nb_genotypes_amr, nb_years,nb_genotypes_vstat] ; // pred numbers
	real pred_absolute_freq_amr[nb_countries,nb_genotypes_amr, nb_years,nb_genotypes_vstat] ; // pred numbers
	
	// real pred_number_non_ref_amr[nb_countries,nb_genotypes_amr-1, nb_years, nb_genotypes_vstat] ; // pred numbers non ref
	// real pred_number_ref_amr[nb_countries, nb_years,nb_genotypes_vstat] ; // pred numbers ref
	
	simplex[nb_genotypes_amr] f0_amr[nb_countries,nb_genotypes_vstat];
  
	
	//////////////////VACCINE TYPE PARAMETERS
	real<lower=-5, upper=5> fitness_genotypes_vector_pre_vacc_vstat[nb_genotypes_vstat-1, 1];
	real<lower=-5, upper=5> fitness_genotypes_vector_post_vacc_vstat[nb_genotypes_vstat-1, 1];
	
	real<lower=0, upper=1> pred_paired_freq_vstat[nb_countries,nb_genotypes_vstat-1, nb_years] ; // pred frequencies
	real pred_paired_freq_tmp_vstat[nb_countries,nb_genotypes_vstat, nb_years] ; // pred numbers
	real pred_absolute_freq_vstat[nb_countries,nb_genotypes_vstat, nb_years] ; // pred numbers
	
	// real pred_number_non_ref_vstat[nb_countries,nb_genotypes_vstat-1, nb_years] ; // pred numbers non ref
	// real pred_number_ref_vstat[nb_countries, nb_years] ; // pred numbers ref
	
	simplex[nb_genotypes_vstat] f0_vstat[nb_countries];

	real pred_absolute_freq_total[nb_countries,nb_genotypes_total, nb_years] ; // pred numbers
	real pred_number_total[nb_countries,nb_genotypes_total, nb_years] ; // pred numbers ref
  // real<lower=0, upper=10> overdispersion_total;


	
//calculating starting frequency from the total
//   for(k in 1:nb_countries){
// 	  for(j in 1:nb_genotypes_amr){
// 	  f0_amr[k,j] = f0_total[k,j]+f0_total[k,j+2];
// 	  //adjusting proportions
// // 	  if(j<nb_genotypes_amr){
// //     f0_amr[k,j] = ((f0_total[k,j]+f0_total[k,j+2])*(freq_per_class[j]/freq_per_class_all));
// //     }
// //     if(j==nb_genotypes_amr){
// //     f0_amr[k,j] = ((f0_total[k,j]+f0_total[k,j+2])*(freq_per_class[j]/freq_per_class_all));
// //     }
//       
// 	  }
//   }
  
// Nvt-s, nvt-r, vt-s, vt-r
///calculating starting frequency from the total with nvt and vt
for(k in 1:nb_countries){
   //computing for nonvaccine types
  for(j in 1:nb_genotypes_amr){
    f0_amr[k,1,j] = f0_total[k,j]/sum(f0_total[k,1:2]);
  }
   //computing for vaccine types
  for(j in 1:nb_genotypes_amr){
    f0_amr[k,2,j] = f0_total[k,j+2]/sum(f0_total[k,3:4]);
  }
}
	
	//calculate starting frequency vstat from total
	for(k in 1:nb_countries){
	  for(j in 1:nb_genotypes_vstat){
	    if(j<nb_genotypes_vstat){
	f0_vstat[k,j] = f0_total[k,j]+f0_total[k,j+1];
	  }
	if(j==nb_genotypes_vstat){
	  f0_vstat[k,j] = f0_total[k,j+1]+f0_total[k,j+2];
	    }
	  }
	}
	
	// /adjusting proportions
	// for(k in 1:nb_countries){
	//   for(j in 1:nb_genotypes_vstat){
	//     if(j<nb_genotypes_vstat){
	// f0_vstat[k,j] = ((f0_total[k,j]+f0_total[k,j+1])*(freq_per_class[j+2]/freq_per_class_all));
	//   }
	// if(j==nb_genotypes_vstat){
	//   f0_vstat[k,j] = ((f0_total[k,j+1]+f0_total[k,j+2])*(freq_per_class[j+2]/freq_per_class_all));
	//     }
	//   }
	// }
	
		////////////////////////AMR
	
	// construct fitness vector
	fitness_genotypes_vector_pre_vacc_amr = fitness_genotypes_pre_vacc_amr;
	fitness_genotypes_vector_post_vacc_amr = fitness_genotypes_post_vacc_amr;
  
  // compute predicted numbers
// 	for (k in 1:nb_countries){
// 	  // compute pred fequencies
// 	  for (j in 1:(nb_genotypes_amr-1)){///this is estimating a single parameter for amr
// 	    pred_paired_freq_amr[k,j,] = logistic_model(f0_amr[k,j]/(f0_amr[k,j]+f0_amr[k,nb_genotypes_amr]), fitness_genotypes_vector_pre_vacc_amr[j,1], fitness_genotypes_vector_post_vacc_amr[j,1], nb_years, yearF0[k], yearIntroduction[k]);
// 	  }
//     for(l in 1:nb_years){
//  	  for (j in 1:(nb_genotypes_amr-1)){///this is estimating a single parameter for amr
//             pred_paired_freq_tmp_amr[k,j,l] = (pred_paired_freq_amr[k,j,l])*10000/(1-pred_paired_freq_amr[k,j,l]);
//             if(pred_paired_freq_tmp_amr[k,j,l] < 0) pred_paired_freq_tmp_amr[k,j,l] = 0;
//     	  }
//     	  pred_paired_freq_tmp_amr[k,nb_genotypes_amr,l] = 10000;
//     	  for(j in 1:nb_genotypes_amr){
//     	    pred_absolute_freq_amr[k,j,l] = 0;
//     	    pred_absolute_freq_amr[k,j,l] = pred_paired_freq_tmp_amr[k,j,l]/sum(pred_paired_freq_tmp_amr[k,,l]);
//     	    if(pred_absolute_freq_amr[k,j,l]<0) pred_absolute_freq_amr[k,j,l] = 0;
//     	  }
//     }
// 	}
	
	
	  // compute predicted numbers
	for (k in 1:nb_countries){
	  // compute pred fequencies
	  for (v in 1:nb_genotypes_vstat){ 
	  for (j in 1:(nb_genotypes_amr-1)){///this is estimating a single parameter for amr
	  // for (j in 1:(nb_genotypes_amr)){//this is allowing a two parameters for amr
	  	pred_paired_freq_amr[k,j,,v] = logistic_model(f0_amr[k,v,j]/(f0_amr[k,v,j]+f0_amr[k,v,nb_genotypes_amr]), fitness_genotypes_vector_pre_vacc_amr[v,1], fitness_genotypes_vector_post_vacc_amr[v,1], nb_years, yearF0[k], yearIntroduction[k]);
	    // pred_paired_freq_amr[k,j,,v] = logistic_model(f0_amr[k,j,v]/(f0_amr[k,j,v]+f0_amr[k,nb_genotypes_amr,v]), fitness_genotypes_vector_pre_vacc_amr[v,1], fitness_genotypes_vector_post_vacc_amr[v,1], nb_years, yearF0[k], yearIntroduction[k]);
	  }
    for(l in 1:nb_years){
  	  for (j in 1:(nb_genotypes_amr-1)){///this is estimating a single parameter for amr
	          // for (j in 1:(nb_genotypes_amr)){//this is allowing a two parameters for amr
            pred_paired_freq_tmp_amr[k,j,l,v] = (pred_paired_freq_amr[k,j,l,v])*10000/(1-pred_paired_freq_amr[k,j,l,v]);
            if(pred_paired_freq_tmp_amr[k,j,l,v] < 0) pred_paired_freq_tmp_amr[k,j,l,v] = 0;
    	  }
    	  pred_paired_freq_tmp_amr[k,nb_genotypes_amr,l,v] = 10000;
    	  for(j in 1:nb_genotypes_amr){
    	    pred_absolute_freq_amr[k,j,l,v] = 0;
    	    pred_absolute_freq_amr[k,j,l,v] = pred_paired_freq_tmp_amr[k,j,l,v]/sum(pred_paired_freq_tmp_amr[k,,l,v]);
    	    if(pred_absolute_freq_amr[k,j,l,v]<0) pred_absolute_freq_amr[k,j,l,v] = 0;
    	  }
    }
	}
	}
	
	////////////////////////Vaccine Status
	
	// construct fitness vector
	fitness_genotypes_vector_pre_vacc_vstat = fitness_genotypes_pre_vacc_vstat;
	fitness_genotypes_vector_post_vacc_vstat = fitness_genotypes_post_vacc_vstat;
  
  // compute predicted numbers
	for (k in 1:nb_countries){
	  // compute pred fequencies
	  for (j in 1:(nb_genotypes_vstat-1)){
	    pred_paired_freq_vstat[k,j,] = logistic_model(f0_vstat[k,j]/(f0_vstat[k,j]+f0_vstat[k,nb_genotypes_vstat]), fitness_genotypes_vector_pre_vacc_vstat[j,1], fitness_genotypes_vector_post_vacc_vstat[j,1], nb_years, yearF0[k], yearIntroduction[k]);
	  }
    for(l in 1:nb_years){
  	    for(j in 1:(nb_genotypes_vstat-1)){
  	       // for(j in 2:(nb_genotypes_vstat)){
            // pred_paired_freq_tmp_vstat[k,j,l] = (pred_paired_freq_vstat[k,j,l])*10000/(1-pred_paired_freq_vstat[k,j,l]);
            pred_paired_freq_tmp_vstat[k,2,l] = (pred_paired_freq_vstat[k,j,l])*10000/(1-pred_paired_freq_vstat[k,j,l]);
            // if(pred_paired_freq_tmp_vstat[k,j,l] < 0) pred_paired_freq_tmp_vstat[k,j,l] = 0;
            if(pred_paired_freq_tmp_vstat[k,2,l] < 0) pred_paired_freq_tmp_vstat[k,2,l] = 0;
    	  }
    	  pred_paired_freq_tmp_vstat[k,nb_genotypes_vstat-1,l] = 10000;///Change it so that I compute for the non-ref (VT) first
    	  // pred_paired_freq_tmp_vstat[k,nb_genotypes_vstat,l] = 10000;
    	  for(j in 1:nb_genotypes_vstat){
    	    pred_absolute_freq_vstat[k,j,l] = 0;
    	    pred_absolute_freq_vstat[k,j,l] = pred_paired_freq_tmp_vstat[k,j,l]/sum(pred_paired_freq_tmp_vstat[k,,l]);
    	    if(pred_absolute_freq_vstat[k,j,l]<0) pred_absolute_freq_vstat[k,j,l] = 0;
    	  }
    }
	}

		////////////////////////Vaccine Status and AMR

	// overdispersion
	// overdispersion_total = 1/sqrt(overdispersion_inv_total);

  
  // compute predicted numbers
	for (k in 1:nb_countries){
    for(l in 1:nb_years){
    	    for(v in 1:nb_genotypes_vstat){
    	      for(m in 1:nb_genotypes_amr){
    	    pred_absolute_freq_total[k,((v-1)*(nb_genotypes_amr)+m),l] = 0;
    	    // pred_absolute_freq_total[k,((v-1)*(nb_genotypes_amr)+m),l] = pred_absolute_freq_amr[k,m,l]*pred_absolute_freq_vstat[k,v,l];
    	   pred_absolute_freq_total[k,((v-1)*(nb_genotypes_amr)+m),l] = pred_absolute_freq_amr[k,m,l,v]*pred_absolute_freq_vstat[k,v,l];
    	    if(pred_absolute_freq_total[k,((v-1)*(nb_genotypes_amr)+m),l]<0) pred_absolute_freq_total[k,((v-1)*(nb_genotypes_amr)+m),l] = 0;
    	      }
    	    }

    	  for(j in 1:(nb_genotypes_total)){
    	    pred_number_total[k,j,l] = pred_absolute_freq_total[k,j,l]*(data_total_number_total[l,k]);
    	    if (pred_number_total[k,j,l] <1E-30) {
    		  	pred_number_total[k,j,l] = 1E-30;
    		  }
	      }

    }
	}
}

///////MAKE THE NAMES MATCH
model {
  
  //
   // overdispersion_inv_total ~ normal(0, 1);
  
  
  ///swap index to test
    //estimate fitness for first parameter for amr
    for (i in 1:(nb_genotypes_amr-1)){
    fitness_genotypes_pre_vacc_amr[1,i] ~ cauchy(0, 0.15);
    fitness_genotypes_post_vacc_amr[1,i] ~ cauchy(fitness_genotypes_pre_vacc_amr[1,i], 0.15);
    // fitness_genotypes_post_vacc_amr[1,i] ~ cauchy(0, 0.15);

  }
  //estimate fitness for second parameter for amr
   for (i in 1:(nb_genotypes_amr-1)){
    fitness_genotypes_pre_vacc_amr[2,i] ~ cauchy(0, 0.15);
    fitness_genotypes_post_vacc_amr[2,i] ~ cauchy(fitness_genotypes_pre_vacc_amr[2,i], 0.15);
    // fitness_genotypes_post_vacc_amr[2,i] ~ cauchy(0, 0.15);

  }
  
  //estimate fitness for first parameter for amr
  // for (i in 1:(nb_genotypes_amr-1)){
  //   fitness_genotypes_pre_vacc_amr[i,1] ~ cauchy(0, 0.15);
  //   fitness_genotypes_post_vacc_amr[i,1] ~ cauchy(fitness_genotypes_post_vacc_amr[i,1], 0.15);
  // }
  // //estimate fitness for second parameter for amr
  //  for (i in 1:(nb_genotypes_amr-1)){
  //   fitness_genotypes_pre_vacc_amr[i,2] ~ cauchy(0, 0.15);
  //   fitness_genotypes_post_vacc_amr[i,2] ~ cauchy(fitness_genotypes_post_vacc_amr[i,2], 0.15);
  // }
  
  ///estimate fitness for nvt and vt
    for (i in 1:(nb_genotypes_vstat-1)){
    fitness_genotypes_pre_vacc_vstat[i,1] ~ cauchy(0, 0.15);
    // fitness_genotypes_post_vacc_vstat[i,1] ~ cauchy(fitness_genotypes_pre_vacc_vstat[i,1], 0.15);
    fitness_genotypes_post_vacc_vstat[i,1] ~ cauchy(0, 0.15);
  }
  //
	for (k in 1:nb_countries){
    for(l in 1:nb_years){
  	    for (j in 1:(nb_genotypes_total)){
  	       // for (j in 1:(nb_genotypes_total-1)){
  	      //Vaccine Status and AMR
  	      if(non_zero_country_year_total[l,k]) {
  	 		    // data_genotype_total[j,l,k] ~ neg_binomial_2(pred_number_total[k,j,l], pow(pred_number_total[k,j,l], overdispersion_total));
  	 		    data_genotype_total[j,l,k] ~ poisson(pred_number_total[k,j,l]);
  	      }
	  	  }
      }
	  }
}

generated quantities {
  vector[nb_countries*(nb_genotypes_total)*nb_years-(number_zeros_country_year_total*nb_genotypes_total)] log_lik;
  
  // vector[nb_countries*(nb_genotypes_total)*nb_years-number_zeros_country_year_genotype_total] log_lik;
  // vector[nb_countries*(nb_genotypes_total-1)*nb_years-number_zeros_country_year_genotype_total] log_lik;
  int n;
  n = 1;
  for (k in 1:nb_countries){
	  for(l in 1:nb_years){
	      for (j in 1:(nb_genotypes_total)){
	       // for (j in 1:(nb_genotypes_total-1)){
	        if(non_zero_country_year_total[l,k]) {
    	      log_lik[n] = 1E-30;
    	      // log_lik[n] = neg_binomial_2_lpmf(data_genotype_total[j,l,k]| pred_number_total[k,j,l], pow(pred_number_total[k,j,l], overdispersion_total));
     	      log_lik[n] = poisson_lpmf(data_genotype_total[j,l,k]| pred_number_total[k,j,l]);
  	 		    n=n+1;
	        }
	      }
  	}
	}
}
