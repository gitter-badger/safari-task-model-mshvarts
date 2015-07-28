library(data.table)
library(plyr)

setwd('~/OneDrive/repos/safari-data/')

# fread() does not play nicely with na.strings apparently
tours <- data.table(read.csv('tours.csv', na.strings=c("NaN", "NA")))

# exclude subject 2 for weirdness as per email from Stephanie 7/28
tours <- tours[subject!=2]

# compute RT
tours[,rt:=t_response-t_question]

# trials right now index within sector/session. we want a global one by subj
tours[,trialWithinSector:=trial]
tours[,trial:=1:.N,by="subject"]

genProb <- matrix(c(0.1141, 0.2617, 0.2571, 0.4563, 0.1748, 0.0778, 0.1441, 0.2774, 0.0735, 0.0677, 0.1222, 0.0166, 0.3690, 0.2181, 0.2522, 0.0271, 0.2686, 0.3747, 0.2244, 0.2226), ncol=4, byrow=T)
colSums(genProb) # sum to 1, sanity check

nSectors <- ncol(genProb)
nAnimals <- nrow(genProb)
nSubjects <- length(unique(tours$subject))
nTrials <- max(tours$trial)

# I should have a posterior that is nsubjects x ntrials x nSectors x nAnimals 
getCounts <- function(d){
  counts <- array(dim=c(nTrials+1, nSectors, nAnimals))
  counts[1,,] <- 1 # proper uniform prior
  for (i in 1:nrow(d)){ 
    counts[i+1,,] <- counts[i,,] # copy last time point
    counts[i+1,d[i,]$sector, d[i,]$animal] <- counts[i,d[i,]$sector, d[i,]$animal] + 1 # increment what we saw
  }
  return(counts)
}

# get the posteriors (assuming dirichlet prior/post) for each timepoint
counts <- daply(tours, .(subject), getCounts, .progress="text")
counts <- adply(counts, 1:3, .progress="text") # get this back into a table (one row per sector per subject per trial)

colnames(counts)[2:3] <- c("trial", "sector")

# slightly annoying: adply turns subject, trial and sector into factors. we want them as integers, as per ?factor: 
counts$subject <- as.numeric(levels(counts$subject))[counts$subject] 
counts$trial <- as.numeric(levels(counts$trial))[counts$trial] 
counts$sector <- as.numeric(levels(counts$sector))[counts$sector] 

counts <- melt(counts, id.vars=c("subject","trial","sector"), variable.name="animal", value.name="countSoFar")

# at a first stab, let's assume there's nothing interesting going on a sector other than the one we're in now. so let's just look up the counts into the data
head(tours)
tours <- merge.data.frame(tours, counts, by.x=c("subject", "trial", "sector","question_1"), by.y=c("subject","trial","sector","animal"), sort=F)
colnames(tours)[15] <- "animal1Count"
tours <- merge.data.frame(tours, counts, by.x=c("subject", "trial", "sector","question_2"), by.y=c("subject","trial","sector","animal"), sort=F)
colnames(tours)[16] <- "animal2Count"
head(tours)

tours <- as.data.table(tours)

tours[,trialsThisSector:=1:.N,by="subject,sector"]
tours[,trialsThisSector:=1:.N,by="subject,sector"]

tours[,animal1Prior:=animal1Count/(trialsThisSector+1)]
tours[,animal2Prior:=animal2Count/(trialsThisSector+1)]

tours[,animalLogOdds:=log(animal1Prior/animal2Prior)]

write.csv(tours, "tours_preprocessed.csv", row.names=F)
