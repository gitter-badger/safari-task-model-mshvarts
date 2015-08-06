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
  // global
  int<lower=1> N_SUBJ; 
  int<lower=1> N_SECTORS; 
  int<lower=1> N_ANIMALS; 
  
  // tours
  int<lower=1> TOURS_N; // n observations
  real tours_rt[TOURS_N]; 
  int tours_animalChoiceMade[TOURS_N];
  int tours_sectors[TOURS_N];
  int tours_subjects[TOURS_N];
  int tours_animal1Choice[TOURS_N];
  int tours_animal2Choice[TOURS_N];
  int tours_animalShown[TOURS_N];

  // trials 
  int<lower=1> TRIALS_N; // n observations for trials
  int<lower=1> N_TEST_TRIALS; // n observations for trials
  int trials_animalDistr[TRIALS_N,N_ANIMALS];
  int trials_subjects[TRIALS_N]; 
  int trials_secIDs[TRIALS_N];
  int trials_sector1choice[TRIALS_N];
  int trials_sector2choice[TRIALS_N];

}

transformed data{
    int N_SUBJ_PARAMS;
    real rt_conditional[TOURS_N]; 

    N_SUBJ_PARAMS <- 3; 
    for ( i in 1:TOURS_N){
     if(tours_animalChoiceMade[i] == 1){
      rt_conditional[i] <- tours_rt[i];
    } 
    else {
      rt_conditional[i] <- -tours_rt[i];
    }
  }
}

parameters{
  real softmaxAdjust[N_SUBJ]; 
  real<lower=0> posteriorNoise[N_SUBJ]; 
  real<lower=0> subject_t0[N_SUBJ]; 
  real<lower=0> subject_z[N_SUBJ]; 
  real subject_x0[N_SUBJ]; 
  simplex[N_ANIMALS] subject_tours_posteriors[TOURS_N, N_SUBJ, N_SECTORS]; 
}

model{
  
  real logCountDifference;  // real because counts are noisified
  // for each subject we have their counts
  vector[N_ANIMALS] counts[N_SUBJ, N_SECTORS]; 

  # priors 
  softmaxAdjust ~ normal(1, 10); # a priori assumption that subjects use the counts directly
  posteriorNoise ~ cauchy(0, 2.5); # half-cauchy scale prior concentrates mass on 0, with fat tails

  // initialize the counts  
  for (i in 1:N_SUBJ){
    for (j in 1:N_SECTORS){
      counts[i,j] <- rep_vector(1, N_ANIMALS);
    }
  }
  
  for ( i in 1:TOURS_N ) {
    subject_tours_posteriors[i, tours_subjects[i], tours_sectors[i]] ~ normal(counts[tours_subjects[i], tours_sectors[i]], posteriorNoise[tours_subjects[i]]); 
    # get the log noisified ratio
    logCountDifference <- inv_logit(subject_tours_posteriors[i, tours_subjects[i]][tours_sectors[i]][tours_animal1Choice[i]]-subject_tours_posteriors[i, tours_subjects[i]][tours_sectors[i]][tours_animal2Choice[i]]); 
    # draw the choice based on multiplicatively perturbed difference
    rt_conditional[i] ~ wfpt_bogacz(subject_t0[tours_subjects[i]], logCountDifference, subject_z[tours_subjects[i]], subject_x0[tours_subjects[i]]); 
    # then increment the count
    counts[tours_subjects[i]][tours_sectors[i]][tours_animalShown[i]] <- counts[tours_subjects[i]][tours_sectors[i]][tours_animalShown[i]] + 1; 
  }
  

}

generated quantities {

  real trialChoiceDiff; 
  int predicted_choices[TRIALS_N];
  real subject_trials_posteriors[N_SUBJ, N_SECTORS];
  simplex[N_ANIMALS] final_training_posterior[N_SUBJ, N_SECTORS]; 
  real count_sums[N_SUBJ, N_SECTORS]; 
  // normalize the noisified counts to get the final simplex for each subject/sector
  // initialize norms
  for (i in 1:N_SUBJ){
    for (j in 1:N_SECTORS){
      count_sums[i,j] <- 0; 
    }
  }
  # compute sums
  for (i in 1:N_SUBJ){
    for (j in 1:N_SECTORS){
      for (k in 1:N_ANIMALS){
        count_sums[i, j] <- count_sums[i,j] + subject_tours_posteriors[TOURS_N,i,j,k]; 
      }
    }
  }  
  # normalize 
  for (i in 1:N_SUBJ){
    for (j in 1:N_SECTORS){
      final_training_posterior[i, j] <- subject_tours_posteriors[TOURS_N,i,j] / count_sums[i,j]; 
    }
  }
  
  # now we can compute the posterior given our prior (given by the ormalized counts and )
  for (i in 1:TRIALS_N){
    # compute posteriors
    for (j in 1:N_SECTORS){
      subject_trials_posteriors[trials_subjects[i], j] <- multinomial_log(trials_animalDistr[i], final_training_posterior[trials_subjects[i],j]);
    }    
    // look up into posterior based on choices made
    trialChoiceDiff <- subject_trials_posteriors[trials_subjects[i], trials_sector1choice[i]] - subject_trials_posteriors[trials_subjects[i], trials_sector2choice[i]];
    predicted_choices[i] <- bernoulli_rng(inv_logit(softmaxAdjust[trials_subjects[i]] * trialChoiceDiff));
  }
}
