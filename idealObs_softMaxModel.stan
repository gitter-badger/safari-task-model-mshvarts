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


parameters{
  real softmaxAdjustHyperPrior; 
  real softmaxAdjust[N_SUBJ]; 
  real<lower=0> posteriorNoise[N_SUBJ]; 
  simplex[N_ANIMALS] subject_tours_posteriors[TOURS_N, N_SUBJ, N_SECTORS]; 
}

model{
  
  real countDifference;  // real because counts are noisified
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
    # the posterior is a noisified count (this is the simplest assumption that lets us treat the posterior as a random variable)
    subject_tours_posteriors[i, tours_subjects[i], tours_sectors[i]] ~ normal(counts[tours_subjects[i], tours_sectors[i]], posteriorNoise[tours_subjects[i]]); 
    # get the log noisified ratio
    countDifference <- subject_tours_posteriors[i, tours_subjects[i]][tours_sectors[i]][tours_animal1Choice[i]]-subject_tours_posteriors[i, tours_subjects[i]][tours_sectors[i]][tours_animal2Choice[i]]; 
    # draw the choice based on multiplicatively perturbed difference
    tours_animalChoiceMade[i] ~ bernoulli_logit(softmaxAdjust[tours_subjects[i]]*countDifference); 
    # then increment the count
    counts[tours_subjects[i]][tours_sectors[i]][tours_animalShown[i]] <- counts[tours_subjects[i]][tours_sectors[i]][tours_animalShown[i]] + 1; 
  }
  

}

generated quantities {

  real trialChoiceDiff; 
  int predicted_choices[TRIALS_N];
  real subject_trials_posteriors[N_SUBJ, N_SECTORS];
  simplex[N_ANIMALS] final_training_posterior[N_SECTORS, N_SECTORS]; 
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
