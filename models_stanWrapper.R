library(rstan)
library(data.table)


tours <- fread('./tours_preprocessed.csv')

# drop timeouts
tours <- tours[!is.na(response)]
tours[,rt:=t_response-t_question]
tours[,response01:=ifelse(response==1, 1, 0)]

trials <- fread('./trials_preprocessed.csv')

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

stanmodel <- stan(file="idealObs_softMaxModel.stan", data=standata, chains=1, iter=1)

fit_idealObs_softmax <- stan(fit=stanmodel, data=standata, chains=3, iter=300, warmup=100, cores=3)