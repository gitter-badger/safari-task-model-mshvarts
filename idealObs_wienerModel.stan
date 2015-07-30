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
     	// if rt < t0, return something really small (otherwise degeneracy ensues)
     	if (rt_pos < t0){
     		return log(1e-10);
     	}
     	// transform params, flip based on conditional rt
     	if (rt>0) { 
     		beta <- x0;
     		delta <- a; 
     	} 
     	else {
     		beta <- 1 - x0;
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
       # print("dwiener(",rt_pos,",",alpha,",",t0,",",beta,",",delta,")");
       return(wiener_log(rt_pos, alpha, t0, beta, delta)); 
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
    int N_SUBJ_PARAMS;
    real rt_conditional[N]; 

    N_SUBJ_PARAMS <- 3; 
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
  vector[N_SUBJ_PARAMS] subjGrandMeans; 
  vector<lower=0> [N_SUBJ_PARAMS] subject_coefs[N_SUBJ]; // subject coefs in a vector so we can define covariance mat
  cholesky_factor_corr[N_SUBJ_PARAMS] L_Omega_subj; // cholesky factorization of subject effects correlation matrix
  vector<lower=0>[N_SUBJ_PARAMS] tau_subj; // scale vector for subject effects correlation matrix

}

transformed parameters{ 
  real subject_x0[N_SUBJ];
  real<lower=0> subject_z[N_SUBJ]; 
  real<lower=0> subject_t0[N_SUBJ];
  real x0_grandMean; 
  real z_grandMean; 
  real t0_grandMean; 

  for (i in 1:N_SUBJ){
    subject_x0[i] <- subject_coefs[i][1];
    subject_z[i] <- subject_coefs[i][2];
    subject_t0[i] <- subject_coefs[i][3];
  }

  x0_grandMean <- subjGrandMeans[1]; 
  z_grandMean <- subjGrandMeans[2]; 
  t0_grandMean <- subjGrandMeans[3]; 

}
model{
  // declarations 
  matrix[N_SUBJ_PARAMS,N_SUBJ_PARAMS] Sigma_subj; // covariance matrix for subject effects, computed from corr matrix and scale
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
  // subject ranefs
  L_Omega_subj ~ lkj_corr_cholesky(2); // prior on the cholesky factorization of the correlation matrix of the subject ranefs
  tau_subj ~ cauchy(0, 2.5); // prior on scale on the correlation matrix of the subject ranefs
  Sigma_subj <- diag_pre_multiply(tau_subj, L_Omega_subj); // diag(tau) * Omega, as per manual p59
  subject_coefs ~ multi_normal_cholesky(subjGrandMeans, Sigma_subj); 

  for ( i in 1:N ) {
   # get the log count ratio (means no need to normalize!)
   logCountRatio <- log(counts[subjects[i]][sectors[i]][animal1Choice[i]]/counts[subjects[i]][sectors[i]][animal2Choice[i]]); 
   rt_conditional[i] ~ wfpt_bogacz(subject_t0[subjects[i]], logCountRatio, subject_z[subjects[i]], subject_x0[subjects[i]]); 
   counts[subjects[i]][sectors[i]][animalShown[i]] <- counts[subjects[i]][sectors[i]][animalShown[i]] + 1; 
 }

}
