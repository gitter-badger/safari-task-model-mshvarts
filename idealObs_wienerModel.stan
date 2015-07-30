functions {
	    /**
     * Compute WFPT density using the parameterization of 
     * Bogacz et al. 2006
     * @param rt Response time. 
     * @param choice Choice identity (0 for lower, 2 for upper)
     * @param t0 Nondecision time. 
     * @param a Drift rate (absolute)
     * @param z Threshold (absolute)
     * @param x0 Startpoint (absolute)
     * @return density of WFPT with parameters t0, a, z, x0 evaluated at conditional RT given by rt, choice
     */
     real wfpt_bogacz_log(real rt, real t0, real a, real z, real x0) {
     	// wiener_log expects: 
     	// alpha: boundary separation (delta between boundaries, aka 2z)
     	// tau: nondecision time
     	// beta bias (0..1, scales between ±z)
     	// delta drift (unconstrained)
     	real rt_pos;
     	real alpha;
     	real tau;
     	real beta; 
     	real delta; 
     	rt_pos <- fabs(rt); // get rid of the weird negative RT convention now
     	alpha <- 2 * z; 
     	tau <- t0;
     	// if rt < t0, return something really small (otherwise degeneracy ensues)
     	if (rt_pos < t0){
     		return log(1e-10);
     	}
     	// transform params, flip based on conditional rt
     	if (rt>0) { 
     		beta <- (x0 + z) / (2 * z);
     		delta <- a; 
     	} 
     	else {
     		beta <- 1 - (x0 + z) / (2 * z);
     		delta <- -a; 
     	}
			// if startpoint is above/below threshold, set very close to threshold
			// (if this is too close, the mass is a delta function on t0, I think .999 is reasonable)
     	if (beta >= 1){
     		beta <- 0.999;
     	} 
     	if (beta <= 0){
     		beta <- 0.0001; 
     	}

     	# print(rt);
     	# print(t0);
     	# print(a);
     	# print(z);
     	# print(x0);
     	# print("alpha=",alpha, "tau=",tau,"beta=",beta,"delta=",delta,"rt_pos=",rt_pos);
     	# print("dwiener(",rt_pos,",",alpha,",",tau,",",beta,",",delta,")");
     	return(wiener_log(rt_pos, alpha, tau, beta, delta)); 
     }

   }

   data{
   	int<lower=1> N; // n observations
   	int<lower=1> N_SUBJ; 
   	int<lower=1> N_SECTORS; 
   	int<lower=1> N_ANIMALS; 
   	real rt[N]; 
   	int animalChoiceMade[N];
   	int sectors[N];
   	int subjects[N];
   	int animal1Choice[N];
   	int animal2Choice[N];
   	int animalShown[N];

   }
   transformed data{
   	real rt_conditional[N]; 
   	for ( i in 1:N){
   		if(animalChoiceMade[i] == 1){
   			rt_conditional[i] <- rt[i];
   		} 
   		else {
   			rt_conditional[i] <- -rt[i];
   		}
   	}
   }

   parameters{
   	real<lower=0> subj_a_mult[N_SUBJ]; // multiplier on the log probability ratio to get drift
   	# real<lower=0> subj_a_baseline[N_SUBJ]; 
   	// real x0_mult[N_SUBJ]; // multiplier on the log probability ratio to get startpoint
   	real<lower=0> subj_x0[N_SUBJ]; // NDT is positive
   	real<lower=0> subj_z[N_SUBJ]; // thresholds are positive (symmetric)
   	real<lower=0> subj_t0[N_SUBJ]; // t0s are positive
   }

   model{
   	// declarations
   	real a; // compute drift from log count ratio and multiplier
   	real logCountRatio; 
   	
   	// for each subject we have their counts (in non ideal obs this should be a parameter)
   	vector[N_ANIMALS] counts[N_SUBJ, N_SECTORS]; 
   	// initialize the counts
   	for (i in 1:N_SUBJ){
   		for (j in 1:N_SECTORS){
   			for (k in 1:N_ANIMALS){
   				counts[i][j][k] <- 1;  // later this can be the prior concentration
   			}
   		}
   	}    	
   	// weak priors
   	subj_a_mult ~ normal(0, 10); 
   	subj_x0 ~ normal(0, 10); 
   	subj_z ~ normal(0, 10); 

   	// pretty stong prior on ndt
   	subj_t0 ~ normal(0.3, 1); 

   	for ( i in 1:N ) {
   		# get the log count ratio (means no need to normalize!)
   		logCountRatio <- log(counts[subjects[i]][sectors[i]][animal1Choice[i]]/counts[subjects[i]][sectors[i]][animal2Choice[i]]); 
   		a <- subj_a_mult[subjects[i]]*logCountRatio; 
   		# x0 <-x0_mult[subjects[i]]*counts[subjects[i]][sectors[i]][animals[i]]; 
   		rt_conditional[i] ~ wfpt_bogacz(subj_t0[subjects[i]], a, subj_z[subjects[i]], subj_x0[subjects[i]]); 
   		counts[subjects[i]][sectors[i]][animalShown[i]] <- counts[subjects[i]][sectors[i]][animalShown[i]] + 1; 
   	}

   }
