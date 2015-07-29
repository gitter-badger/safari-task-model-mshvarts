library(data.table)
library(RWiener)
library(minqa)
library(doMC)

registerDoMC(cores=7) # ncores... 7 for bigmem, fewer for laptop

d <- fread('tours_preprocessed.csv')

# drop timeouts
d <- d[!is.na(response)]

d[,bound:=ifelse(response==1,"upper","lower")]

getWienerLik <- function(params, rt, animalLogOdds, bound){
  alpha <- params[1] # boundary separation
  tau <- params[2] # nondecision time
  beta <- params[3] * animalLogOdds # bias
  # beta is sort of funky: it actually corresponds to a value of ± alpha, but dwiener expects it to be in the range of 0..1 and scales it against alpha internally. Slightly hackish inverse logit here (since we have more than two options so it's not real log odds)
  beta <- exp(beta) / (1+exp(beta))
  delta <- params[4] * animalLogOdds # drift
  negLik <- -dwiener(rt, alpha, tau, beta, delta, resp=bound, give_log=T)
  negLik[is.infinite(negLik)] <- 10000 # very bad
  return(sum(negLik))
}
getWienerLik(c(2, 0.3, 0.5, 1), d$rt, d$animalLogOdds, d$bound)

fitSubject <- function(d,lower = c(0.001, 0.01, -10, -10), upper = c(10, 5, 10, 10), nrestarts = 10){
  cat("Fitting subject %s...", unique(d$subject))
  res <- bobyqa(runif(4, min=lower,max=upper), getWienerLik, lower=lower, upper=upper,control=list(), d$rt, d$animalLogOdds, d$bound)
  bestsofar <- c(res$fval, res$par, res$ierr)
  for (i in 1:nrestarts){
    cat("Restart %s/%s...", i, nrestarts)
    res <- bobyqa(runif(4, min=lower,max=upper), getWienerLik, lower=lower, upper=upper,control=list(), d$rt, d$animalLogOdds, d$bound)
    if (res$fval < bestsofar[1]) bestsofar <- c(res$fval, res$par, res$ierr)
  }
  return(bestsofar)
}

allFits <- ddply(d, .(subject), fitSubject, .parallel=T)
cat("Saving...")
write.csv(allFits, "idealObs_wienerModel_maxLikFits.csv")