# Format Event Log

#Load libraries
library(stringr)
library(lubridate)

#Load data
file = file.choose()
eventlog=read.csv(file)
date = '2021-08-19'

#Format time column
eventlog$Time.stamp.no.ms=strtrim(eventlog$Time.stamp,8)

#Save
dir <- dirname(file)
write.csv(eventlog, file=paste(dir, '/EVENTLOG_', date,'_formatted.csv',sep=""), row.names = F)
 