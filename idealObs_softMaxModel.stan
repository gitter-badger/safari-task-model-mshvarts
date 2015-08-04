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
  real subjGrandMean; 
  real<lower=0> subjHyperVariance;
  real<lower=0> softmaxAdjust[N_SUBJ]; 
  real posteriorNoise[N_SUBJ]; 
  simplex[N_ANIMALS] subject_tours_posteriors[N_SUBJ, N_SECTORS]; 

}

model{
  real countDifference;  // real to account for extension to non-count posteriors later
  // for each subject we have their counts (in non ideal obs this should be a parameter)
  vector[N_ANIMALS] counts[N_SUBJ, N_SECTORS]; 
  real count_sums[N_SUBJ, N_SECTORS]; 
  // initialize the counts  
  for (i in 1:N_SUBJ){
   for (j in 1:N_SECTORS){
    counts[i,j] <- rep_vector(1, N_ANIMALS);
  }
}

for ( i in 1:TOURS_N ) {
 # get the log count ratio (means no need to normalize!)
 countDifference <- counts[tours_subjects[i]][tours_sectors[i]][tours_animal1Choice[i]]-counts[tours_subjects[i]][tours_sectors[i]][tours_animal2Choice[i]]; 
 tours_animalChoiceMade ~ bernoulli_logit(softmaxAdjust[tours_subjects[i]]*countDifference); 
 counts[tours_subjects[i]][tours_sectors[i]][tours_animalShown[i]] <- counts[tours_subjects[i]][tours_sectors[i]][tours_animalShown[i]] + 1; 
}
# now compute what our final posterior is (normalized counts)
// initialize norms
for (i in 1:N_SUBJ){
 for (j in 1:N_SECTORS){
  counts[i,j] <- rep_vector(0, N_ANIMALS);
  count_sums[i,j] <- 0; 
}
}
# compute sums
for (i in 1:N_SUBJ){
  for (j in 1:N_SECTORS){
    for (k in 1:N_ANIMALS){
      count_sums[i, j] <- count_sums[i,j] + counts[i,j,k]; 
    }
  }
}  
# normalize 
for (i in 1:N_SUBJ){
  for (j in 1:N_SECTORS){
      counts[i, j] <- counts[i,j] / count_sums[i,j]; 
    }
}
# now sample from our end-of-training posterior
for (i in 1:N_SUBJ){
  for (j in 1:N_SECTORS){
    subject_tours_posteriors[i,j] ~ normal(counts[i,j], posteriorNoise[i]); 
  }
}
}

generated quantities {

  real trialChoiceDiff; 
  int predicted_choices[TRIALS_N];
  real subject_trials_posteriors[N_SUBJ, N_SECTORS];
  
  
  
  # now we can compute the posterior given our prior (given by the ormalized counts and )
  for (i in 1:TRIALS_N){
    # compute posteriors
    subject_trials_posteriors[trials_subjects[i], trials_secIDs[i]] <- categorical_log(trials_animalDistr[i], subject_tours_posteriors[trials_subjects[i],trials_secIDs[i]]);
    // look up into posterior based on choices made
    trialChoiceDiff <- subject_trials_posteriors[trials_subjects[i], trials_sector1choice[i]] - subject_trials_posteriors[trials_subjects[i], trials_sector2choice[i]];
    predicted_choices[i] <- bernoulli_rng(softmaxAdjust[trials_subjects[i]] * trialChoiceDiff);
  }
}
