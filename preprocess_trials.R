library(rstan)
library(data.table)
library(plyr)

setwd('~/OneDrive/repos/safari-task-model/')
d <- fread('../safari-data/tours.csv')
d.trials <- fread('../safari-data/trials.csv')
# for each test we want the combined evidence (counts) and the response, for each sequence. but sequences are not marked :/

d.trials[,prevTrial:=shift(trial, 1L, type="lag")] 
d.trials[,prevSubject:=shift(subject, 1L, type="lag")] 
d.trials[,isNewSector:=prevTrial>=trial]
d.trials[,isNewSubject:=prevSubject!=subject]

# ugly ugly loop! takes ~5 seconds :(. make me less ugly! 
pb <- txtProgressBar(min=1, max=nrow(d.trials), style=3)
s <- 1
d.trials[1,secId:=1]
for (i in 2:nrow(d.trials)){
  setTxtProgressBar(pb, i)
  if (d.trials[i,isNewSector]) s <- s + 1
  if (d.trials[i,isNewSubject]) s <- 1
  d.trials[i,secId:=s]
}
close(pb)

# for each test we want the combined evidence (counts) and the response, for each sequence. but sequences are not marked :/
trialSummary <- d.trials[,sum(.N),by="secId,subject,animal"]
# we want a column per animal
trialSummary <- dcast(trialSummary, secId+subject~animal, value.var="V1")
colnames(trialSummary)[3:7] <- paste("animal", 1:5, sep="")
nasTo0 <- colnames(trialSummary)[3:7]
# this is fast but obtuse. What it does is select from trialSummary the rows where a column in nasTo0 is NA, set that column for those rows to 0
for (nat0 in nasTo0){
  trialSummary[is.na(get(nat0, trialSummary)),nat0:=0,with=F]
}
head(trialSummary)

trialResponseSummary <- d.trials[trial==1]
setkey(trialResponseSummary, secId, subject)

trialSummary <- trialResponseSummary[trialSummary] # combine

# now trialSummary has one row per trial, just what we want. 
write.csv(trialSummary, 'trials_preprocessed.csv')
