// A Stan model for estimating population growth rates for Indonesia leatherback turtle nests
// at two nesting beaches. 
//
// The state model is normal and the observation model is normal. The probability of nesting in a month
// is modeled with a discrete Fourier series - p[month] sums to 1.0 over each 12 month period.
// To normalize the monthly probability (i.e., they have to sum to 1.0 over 12 months), each
// monthly value has to be divided by the sum.  Constant parts (sum(cos(2*pi*m/P)) and sum(sin(2*pi*m/P),
// where sum is over all m (months) and P is period (6 months for Wermon and 12 months for JM)),
// need to be supplied as data as C_cos and C_sin. Also, C0 is a constant to make all Fourier 
// series values to be > 0.  I think C0 >= 2 should do... (not really... 10 did the trick)

// I used quite informative prior for X[1] ~ dnorm(5, 0.1) although variance of 10 should 
// be wide enough in the log scale.  I used similar priors for the slopes: dnorm(0, 0.1).
// For all standard deviations (process and observation), I used uniform between 0 and 20, 
// which should be wide enough in the log scale. 

// I didn't see any reason to make the priors flatter because that would make convergence difficult
// and did not gain anything in return. 

data {
	int n.years<lower = 1>;
	int n.months<lower = 1>;
	int n.states<lower = 1>;
	int<lower=1> Z[n.states, n.states];  // this may not be needed

	//int<lower = 1> n.obs[n.states];      // numbers of observations
	//int<lower = 1> n.miss[n.states];     // numbers of missing data
	//int<lower = 1, upper = max(n.obs)> idx_obs[n.states, max(n.obs)];
	//int<lower = 1, upper = max(n.miss)> idx_miss[n.states, max(n.miss)];

	real y[n.years, n.months, n.states];   // log(observed # nests per month)
	int m[n.years, n.months];              // months index 1:12 repeats n.years times
	int periods[n.states];        // assumed periodicity of nesting at each location

	real u_mean = 0;
	real u_sd = 1;
}

parameters {

	real<lower=0, upper = 10> sigma.Q[n.states, n.months];       // standard deviation of state process
	real<lower=0, upper = 10> sigma.R[n.states];       // standard deviation of location-specific observation process
	real<lower=0, upper = 10> U[n.states];             // log population growth rates 
	real<lower=0, upper = 10> N[n.states, n.years]     // log total population size
	real p.beta.cos[n.states];
	real p.beta.sin[n.states];
	real C0[n.states];
	real C_cos[n.states];
	real C_sin[n.states];
}

transformed parameters {
	real p.const[n.states, n.months];  // constant part of discrete Fourier series
	real <lower = 0, upper = 1> p[n.states, n.months];  // proportion of nests per month

	for (j in 1:n.states){
        for (k in 1:n.months){
            p.const[j, k] <-  2 * pi * k / period[j]
            p[j, k] <- (C0[j] + p.beta.cos[j] * cos(p.const[j,k]) + p.beta.sin[j] * sin(p.const[j,k]))/(n.months * C0[j] + p.beta.cos[j] * C_cos[j] + p.beta.sin[j] * C_sin[j])
            
        }
    }

	
}

model
{
	// Priors     
	for (j in 1:n.timeseries){
        for (k in 1:n.months){	
			sigma.Q[j,k] ~ dgamma(2, 0.5);
		}
	}

	for (j in 1:n.timeseries){
		sigma.R[j] ~ dgamma(2, 0.5);
		p.beta.cos[j] ~ dnorm(0, 1);
    	p.beta.sin[j] ~ dnorm(0, 1);
    }

    sigma.N ~ dgamma(2, 0.5);

	// ####  Initial states
    for(i in 1:n.states) {
        U[i] ~ dnorm(u_mean, u_sd);
            
        N0[i] ~ dnorm(N0_mean[i], N0_sd[i]);    // prior on init state
        predN[i,1] <- N0[i] + U[i];  
        N[i,1] <- predN[i,1];

    }
   
    // R is the variance of the observation process (Y)
    for(j in 1:n.timeseries) {
        
        for (t in 1:n.months){
           predX[1,t,j] <- log(p[j, m[1,t]]) + N[j,1];
           X[1,t,j] ~ normal(predX[1,t,j], sigma.Q[j,t]);

           // observation
           predY[1,t,j] <- inprod(Z[j, 1:n.states], X[1,t,]);
           
           y[1,t,j] ~  normal(predY[1,t,j], sigma.R[j]);
           
           loglik[1,t,j] <- normal_lpdf(y[1,t,j] | X[1,t,j], sigma.R[j]);   // log likelihood computed here
           
        }       

    }

    //  ####  End of initial states ####

	for (s in 2:n.years){
        for (i in 1:n.states){
            N[i, s] ~ dnorm(pred.N[i, s], sigma.N);
            pred.N[i, s] <- U[i] + N[i, s-1];
        }
            
        for (j in 1:n.timeseries){
            for (t in 1:n.months){
                predX[s,t,j] <- log(p[j, m[s,t]]) + N[j,s]; 
                X[s,t,j] ~ normal(predX[s,t,j], sigma.Q[j,t]);

                // observation
                predY[s,t,j] <- inprod(Z[j, 1:n.states], X[s,t,]);
                y[s,t,j] ~ normal(predY[s,t,j], sigma.R[j]);
                
                loglik[s,t,j] <- normal_lpdf(y[s,t,j] | X[s,t,j], sigma.R[j]);
                
            }       

        }

    }



}

generated quantities {
  // log_lik is for use with the loo package
  vector[Y] log_lik;
  for(i in 1:Y) {
  	log_lik[i] = lognormal_lpdf(H[i] | mu[i], sigma);
  }
}

