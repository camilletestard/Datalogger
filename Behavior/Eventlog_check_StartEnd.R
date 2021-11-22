#Eventlog_check_StartEnd.R
#This script allows to check whether there is an 'end' for every 'start' of behavior in the eventlogs.
# Camille Testard - November 2021

#Load libraries: 
library(timeline)
library(ggplot2) 
library(utils)
library(xlsx)

#Load data:
file = file.choose() # chose the formatted behavior file
log = read.xlsx(file, sheetIndex = 1)

#Get all unique behaviors
behavs = as.character(unique(log$Behavior))
behavs = behavs[!(is.na(behavs) | behavs=="")] #remove blank behaviors

#Only keep behaviors, remove session events
to_remove = cbind("Started recording", "Camera Sync", "Stopped recording")
if (any(is.na(match(to_remove, behavs)))) { #if one of the session event is missing (should only be the case for one session)
  to_remove = to_remove[!is.na(match(to_remove, behavs))]
}
behavs = behavs[-match(to_remove, behavs)] #only keep behaviors

#Split all behavior events by start and end:
new_log=data.frame(); b=1
for (b in 12:length(behavs)){ #For all behaviors
  
  data=log[log$Behavior==behavs[b],c("Time", "Behavior", "Start.end")] #Only keep useful columns
  if (length(which(is.na(data$Behavior)))!=0){data=data[-which(is.na(data$Behavior)),]}# remove rows with empty/NA behavior
  
  interim=split(data,data$Start.end) #Split the data by behavior start and behavior end
  names(interim[[2]])[1]="start.time";names(interim[[1]])[1]="end.time" #rename columns
  
  data = cbind(interim[[1]],interim[[2]]); data=data[,-c(3,5,6)]
  data=data[,c("Behavior","start.time","end.time")]
  data[,c("start.time","end.time")]=data[,c("start.time","end.time")]/1000
  
  new_log=rbind(new_log,data)
}

#Note for Ari: The bit of code above will have an error if there aren't the same
#number of start and end for any one behavior. When an error occurs, check for which behavior it occurs
# and adjust the excel file accordingly (usually an end is missing, or there are 2 starts, use the videos to check). 
# The script stops at the behavior that has an issue, number "b". You can check in the character vector
# 'behavs' which behavior causes issues.
