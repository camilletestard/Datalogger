#Generate behavioral plots with one behavior per line.
library(graphics)
library(matlab)

#Set path
savePath = '~/Dropbox (Penn)/Datalogger/Results/'
dataPath = '~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/'
session_name = 'Amos_2021-07-29'
setwd(paste(dataPath, session_name, sep=""))
labels = read.csv('Labels_per_sec.csv')
behav_categ = names(read.csv('behav_categ.csv'))

block_times = which(diff(labels$block_labels)!=0)

behavs = matrix(0, nrow(labels), length(behav_categ)-1); colnames(behavs) = behav_categ[1:27]
c = 0; empty_behav = vector(); plot.labels = vector()
for (b in 1:length(behav_categ)-1) {
  behavs[, b] = c
  plot.labels[c] = behav_categ[b]
  behavs[labels$behavior_labels_subject == b, b] = c+ 1;
  c = c+2

  if(all(labels$behavior_labels_subject != b)){
    empty_behav = cbind(empty_behav, b)
  }
}

time = replicate(ncol(behavs), 1:nrow(labels))

setwd(paste(savePath, session_name, sep=""))
#png("Behavior_TimeSeries.png")
mar.default <- c(5,10,4,2) + 0.1
matplot(time, behavs, type = "l", axes = F, ylab='')
matlines(ones(c,1)*block_times[1], 1:c, col="black", lwd=1)
matlines(ones(c,1)*block_times[2], 1:c, col="black", lwd=1)
axis(2, at=seq(from=1, to=c-2, by=2)+1, labels=behav_categ[1:27],
     las = 2);
axis(1)
#dev.off()

