// A Stan model for estimating population growth rates for Indonesia leatherback turtle nests
// at two nesting beaches.
//
// The state model is normal and the observation model is normal. The probability of nesting in a month
// is modeled with a discrete Fourier series - p[month] sums to 1.0 over each 12 month period.
// To normalize the monthly probability (i.e., they have to sum to 1.0 over 12 months), each
// monthly value has to be divided by the sum.  Constant parts (sum(cos(2*pi*m/P)) and sum(sin(2*pi*m/P),
// where sum is over all m (months) and P is period (6 months for Wermon and 12 months for JM)),
// need to be supplied as data as C_cos and C_sin. Also, C0 is a constant to make all Fourier 
// series values to be > 0.  I think C0 >= 2 should do (not really but 10 did the trick)

// I used quite informative prior for X[1] ~ normal(5, 0.1) although variance of 10 should 
// be wide enough in the log scale.  I used similar priors for the slopes: normal(0, 0.1).
// For all standard deviations (process and observation), I used uniform between 0 and 20, 
// which should be wide enough in the log scale. 

// I didn't see any reason to make the priors flatter because that would make convergence difficult
// and did not gain anything in return. 

data {
	int<lower = 1> n_years;
	int<lower = 1> n_months;
	int<lower = 1> n_states;
	int<lower=1> Z[n_states, n_states];  // this may not be needed

	//int<lower = 1> n_obs[n_states];      // numbers of observations
	//int<lower = 1> n_miss[n_states];     // numbers of missing data
	//int<lower = 1, upper = max(n_obs)> idx_obs[n_states, max(n_obs)];
	//int<lower = 1, upper = max(n_miss)> idx_miss[n_states, max(n_miss)];

	real y[n_years, n_months, n_states];   // log(observed # nests per month)
	int m[n_years, n_months];              // months index 1:12 repeats n_years times
	int period[n_states];        // assumed periodicity of nesting at each location
	real pi;

}

parameters {

	real<lower=0, upper = 10> sigma_Q[n_states, n_months];       // standard deviation of state process
	real<lower=0, upper = 10> sigma_R[n_states];       // standard deviation of location-specific observation process
	real<lower = 0> sigma_N; 
	real<lower=0> N0[n_states];
	real<lower=0> N0_mean[n_states];
	real<lower=0> N0_sd[n_states];
	real<lower=0, upper = 10> U[n_states];             // log population growth rates 
	real<lower=0, upper = 10> N[n_states, n_years];     // log total population size
	
	real<lower=0, upper = 10> X[n_years, n_states, n_months];     // log total population size
	
	
	real p_beta_cos[n_states];
	real p_beta_sin[n_states];
	real C0[n_states];
	real C_cos[n_states];
	real C_sin[n_states];

}

transformed parameters {
	real p_const[n_states, n_months];  // constant part of discrete Fourier series
	real <lower = 0, upper = 1> p[n_states, n_months];  // proportion of nests per month

	real u_mean = 0;
	real u_sd = 1;

	for (j in 1:n_states){
        for (k in 1:n_months){
            p_const[j, k] =  2 * pi * k / period[j];
            p[j, k] = (C0[j] + p_beta_cos[j] * cos(p_const[j,k]) + p_beta_sin[j] * sin(p_const[j,k]))/(n_months * C0[j] + p_beta_cos[j] * C_cos[j] + p_beta_sin[j] * C_sin[j]);
            
        }
    }

	
}

model
{
	// Priors     
	for (j in 1:n_states){
        for (k in 1:n_months){	
			sigma_Q[j,k] ~ gamma(2, 0.5);
		}
	}

	for (j in 1:n_states){
		sigma_R[j] ~ gamma(2, 0.5);
		p_beta_cos[j] ~ normal(0, 1);
    	p_beta_sin[j] ~ normal(0, 1);
    }

    sigma_N ~ gamma(2, 0.5);

	// ####  Initial states
    for(i in 1:n_states) {
        U[i] ~ normal(u_mean, u_sd);
            
        N0[i] ~ normal(N0_mean[i], N0_sd[i]);    // prior on init state
        //predN[i,1] = N0[i] + U[i];  
        N[i,1] ~ normal(N0[i] + U[i], sigma_N);  //predN[i,1];

    }
   
    // R is the variance of the observation process (Y)
    for(j in 1:n_states) {
        
        for (t in 1:n_months){
           //predX[1,t,j] = log(p[j, m[1,t]]) + N[j,1];
           X[1,t,j] ~ normal(log(p[j, m[1,t]]) + N[j,1], sigma_Q[j,t]);

           // observation
           y[1,t,j] ~  normal(X[1,t,j], sigma_R[j]);
           
        }       

    }

    //  ####  End of initial states ####

	for (s in 2:n_years){
        for (j in 1:n_states){
            //pred_N[j, s] = U[j] + N[j, s-1];
            N[j, s] ~ normal(U[j] + N[j, s-1], sigma_N);
        
            for (t in 1:n_months){
                //predX[s,t,j] = log(p[j, m[s,t]]) + N[j,s]; 
                X[s,t,j] ~ normal(log(p[j, m[s,t]]) + N[j,s], sigma_Q[j,t]);

                // observation
                y[s,t,j] ~ normal(X[s,t,j], sigma_R[j]);                
                
            }       

        }

    }

}

generated quantities {
  // log_lik is for use with the loo package
  real log_lik[n_years, n_states, n_months];
  for (s in 1:n_years){
  	for (j in 1:n_states){
        for (t in 1:n_months){
  			log_lik[s,t,j] = normal_lpdf(y[s,t,j] | X[s,t,j], sigma_R[j]);
  		}
  	}
  }
}

