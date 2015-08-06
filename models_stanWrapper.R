library(rstan)
library(data.table)


tours <- fread('./tours_preprocessed.csv')

# drop timeouts
tours <- tours[!is.na(response)]
tours[,rt:=t_response-t_question]
tours[,response01:=ifelse(response==1, 1, 0)]

trials <- fread('./trials_preprocessed.csv')

# do a subset of subjects
sampled_subjs <- sample(tours$subject, 5)
# sampled_subjs <- unique(tours$subject)
tours <- tours[subject %in% sampled_subjs]
trials <- trials[subject %in% sampled_subjs]

trialsAnimalDistr <- as.matrix(trials[,23:27,with=F])

standata <- list(
  # global
      N_SUBJ = length(unique(tours$subject)),
      N_SECTORS = length(unique(tours$sector)),
      N_ANIMALS = length(unique(tours$animal)),
  # tours
      TOURS_N = nrow(tours),
      tours_rt = tours$rt,
      tours_animalChoiceMade = tours$response01,
      tours_animalShown = tours$animal,
      tours_sectors =  tours$sector,
      tours_subjects =tours$subject,
      tours_animal1Choice = tours$question_1,
      tours_animal2Choice = tours$question_2,
  # trials
      TRIALS_N = nrow(trials),
      N_TEST_TRIALS = length(unique(trials$secId)),
      trials_animalDistr  = trialsAnimalDistr,
      trials_subjects = trials$subject,
      trials_secIDs = trials$secId,
      trials_sector1choice = trials$questions_sector_1,
      trials_sector2choice = trials$questions_sector_2)

# transform subjects back into non-gapped space or we kill stan indexing if we subsampled
standata$tours_subjects <- as.integer(as.factor(standata$tours_subjects))
standata$trials_subjects <- as.integer(as.factor(standata$trials_subjects))

stanmodel <- stan(file="idealObs_softMaxModel.stan", data=standata, chains=1, iter=1)

fit_idealObs_softmax <- stan(fit=stanmodel, data=standata, chains=2, iter=300, warmup=100, cores=2, open_progress=T)