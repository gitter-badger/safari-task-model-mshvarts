library(rstan)
library(data.table)


d <- fread('../safari-data/tours.csv')

# drop timeouts
d <- d[!is.na(response)]
d[,rt:=t_response-t_question]

standata <- list(N = nrow(d),
                 N_SUBJ = length(unique(d$subject)),
                 N_SECTORS = length(unique(d$sector)),
                 N_ANIMALS = length(unique(d$animal)),
                 rt = d$rt,
                 animalChoiceMade = d$response,
                 animalShown = d$animal,
                 sectors = d$sector,
                 subjects = d$subject,
                 animal1Choice = d$question_1,
                 animal2Choice = d$question_2)
                 
  
stanmodel <- stan(file="idealObs_wienerModel.stan", data=standata, chains=1, iter=1)

fit <- stan(fit=stanmodel, data=standata, chains=2, iter=1000, warmup=500)