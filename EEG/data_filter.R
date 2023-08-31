############### PART 0: library, source and pathes ########################

library(eegkit)
library(eegkitdata)
library(dplyr)
library(tidyr)

reading.path <- "D:/VcopularGGM/Codes/Real_data/eeg_full/"
working.path <- "D:/VcopularGGM/Codes/Real_data/"
func.path <- "D:/VcopularGGM/Codes/Functions"
saving.path <- "D:/VcopularGGM/Codes/Real_data/eeg_processed/"

################# PART 1: loading data ######################

########## data from solea ##############
# load(paste(working.path,"X.control.RData",sep=""))
# load(paste(working.path,"X.alcohol.RData",sep=""))

########## read eeg data from eegkitdata package #############

eegS1 <- geteegdata(indir=reading.path, cond="S1", filename="eegfullS1")
save(eegS1, file = 'eeg_raw.rda')

load(file = paste(working.path, 'eeg_raw.rda', sep = ''))
idx.control <- which(eegS1$group=="c")
idx.alcohol <- which(eegS1$group=='a',)

## control group and alcohol group
eegdata.control <- eegS1[idx.control,] %>% distinct()
eegdata.alcohol <- eegS1[idx.alcohol,] %>% distinct()

save(eegdata.control, file = paste(working.path,'eeg_control_raw.rda', sep=''))
save(eegdata.alcohol, file = paste(working.path,'eeg_alcohol_raw.rda', sep=''))

load(file = paste(working.path,'eeg_control_raw.rda', sep=''))
load(file = paste(working.path,'eeg_alcohol_raw.rda', sep=''))

#################### PART 2: processing data ################################

## spread along time
eegdata.control.spread <- spread(eegdata.control, time, voltage)
eegdata.alcohol.spread <- spread(eegdata.alcohol, time, voltage)

## alpha band filtering
time.length <- ncol(eegdata.control.spread) - 5
eegdata.control.alpha <- eegfilter(t(eegdata.control.spread[,-c(1:5)]), Fs=time.length, lower=8, upper=12.5 ) %>%
  t() %>% cbind(eegdata.control.spread[,c(1:5)],.) %>% as.data.frame()
eegdata.alcohol.alpha <- eegfilter(t(eegdata.alcohol.spread[,-c(1:5)]), Fs=time.length, lower=8, upper=12.5 ) %>%
  t() %>% cbind(eegdata.alcohol.spread[,c(1:5)],.) %>% as.data.frame()

## mean over trail
eegdata.control.mean <- aggregate(eegdata.control.alpha[,-c(1:5)], 
                                  list(eegdata.control.alpha$subject, eegdata.control.alpha$channel), mean)
eegdata.alcohol.mean <- aggregate(eegdata.alcohol.alpha[,-c(1:5)], 
                                  list(eegdata.alcohol.alpha$subject, eegdata.alcohol.alpha$channel), mean)

names(eegdata.alcohol.mean)[1:2] <- names(eegdata.control.mean)[1:2] <- c('subject','channel')

subject.list.control <- unique(eegdata.control.mean$subject)
subject.list.alcohol <- unique(eegdata.alcohol.mean$subject)
channel.list <- unique(eegdata.control.mean$channel)

n.control <- length(subject.list.control)
n.alcohol <- length(subject.list.alcohol)
n.channel <- length(channel.list)

eegdata.control.array <- array(data=NA, dim=c(n.control, n.channel, time.length))
eegdata.alcohol.array <- array(data=NA, dim=c(n.alcohol, n.channel, time.length))

for(i in 1:n.control){
  subject.name <- subject.list.control[i]
  for (j in 1:n.channel) {
    channel.name <- channel.list[j]
    eegdata.control.array[i, j, ] <- as.matrix(eegdata.control.mean[eegdata.control.mean$subject==subject.name &
                                                            eegdata.control.mean$channel==channel.name,-c(1:2)])
  }
}

for(i in 1:n.alcohol){
  subject.name <- subject.list.alcohol[i]
  for (j in 1:n.channel) {
    channel.name <- channel.list[j]
    eegdata.alcohol.array[i, j, ] <- as.matrix(eegdata.alcohol.mean[eegdata.alcohol.mean$subject==subject.name &
                                                            eegdata.alcohol.mean$channel==channel.name,-c(1:2)])
  }
}

rownames(eegdata.control.array) <- subject.list.control
rownames(eegdata.alcohol.array) <- subject.list.alcohol
colnames(eegdata.alcohol.array) <- colnames(eegdata.control.array) <- channel.list


save(eegdata.control.array, file = paste(saving.path, 'eeg_control.Rdata', sep=''))
save(eegdata.alcohol.array, file = paste(saving.path, 'eeg_alcohol.Rdata', sep=''))
save(channel.list, file = paste(saving.path, 'nodes_name.Rdata', sep=''))
