#Format_grooming_log.R
#This script re-formats grooming pose and bodypart groomed logs for downstream
#neural data analysis (cross-pose decoding). 

#C. Testard - September 2022
                  
#Load libraries:
library(timeline)
library(ggplot2)
library(utils)
library(xlsx)
library(tidyverse)

#Load data:
mac =1
if (mac == 1) { home = '~'} else
{home = 'C:/Users/GENERAL'}

setwd(paste(home,'/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output', sep=""))
file = file.choose() # chose the formatted behavior file (xls file)
date = as.character(substr(file, nchar(file)-15, nchar(file)-5))
m.monkey = as.character(substr(file, nchar(file)-42, nchar(file)-39))

log <- xlsx::read.xlsx(file, 1)


#Get all unique pose
pose = as.character(unique(log$Pose))
pose = pose[!(is.na(pose) | pose=="")] #remove blank pose

#Get all unique bodyparts
bodypart = as.character(unique(log$Body.part))
bodypart = bodypart[!(is.na(bodypart) | bodypart=="")] #remove blank bodypart

#Get all unique behaviors
behavs = as.character(unique(log$Behavior))
behavs = behavs[!(is.na(behavs) | behavs=="")] #remove blank behaviors

#Only keep behaviors, remove session events
to_remove = cbind("Started recording", "Camera Sync", "Stopped recording")
if (any(is.na(match(to_remove, behavs)))) { #if one of the session event is missing (should only be the case for one session)
  to_remove = to_remove[!is.na(match(to_remove, behavs))]
}
behavs = behavs[-match(to_remove, behavs)] #only keep behaviors

#Split all pose by start and end:
pose_log=data.frame(); p=1
for (p in 1:length(pose)){ #For all behaviors
  
  data=log[log$Pose==pose[p],c("Time", "Pose", "Pose.start.end")] #Only keep useful columns
  if (length(which(is.na(data$Pose)))!=0){data=data[-which(is.na(data$Pose)),]}# remove rows with empty/NA behavior
  
  interim=split(data,data$Pose.start.end) #Split the data by behavior start and behavior end
  names(interim[[2]])[1]="start.time";names(interim[[1]])[1]="end.time" #rename columns
  
  data = cbind(interim[[1]],interim[[2]]); data=data[,-c(3,5,6)]
  data=data[,c("Pose","start.time","end.time")]
  data$start.time=as.numeric(data$start.time); data$end.time=as.numeric(data$end.time)
  data[,c("start.time","end.time")]=data[,c("start.time","end.time")]/1000
  
  pose_log=rbind(pose_log,data)
}

#Split all bodypart by start and end:
bodypart_log=data.frame(); bp=1
for (bp in 1:length(bodypart)){ #For all behaviors
  
  data=log[log$Body.part==bodypart[bp],c("Time", "Body.part", "Pose.start.end")] #Only keep useful columns
  if (length(which(is.na(data$Body.part)))!=0){data=data[-which(is.na(data$Body.part)),]}# remove rows with empty/NA behavior
  
  interim=split(data,data$Pose.start.end) #Split the data by behavior start and behavior end
  names(interim[[2]])[1]="start.time";names(interim[[1]])[1]="end.time" #rename columns
  
  data = cbind(interim[[1]],interim[[2]]); data=data[,-c(3,5,6)]
  data=data[,c("Body.part","start.time","end.time")]
  data$start.time=as.numeric(data$start.time); data$end.time=as.numeric(data$end.time)
  data[,c("start.time","end.time")]=data[,c("start.time","end.time")]/1000
  
  bodypart_log=rbind(bodypart_log,data)
}

#Save to .csv
#output_file = utils::choose.dir(default = "", caption = "Select folder") # choose output directory
# dir <- dirname(output_file)
  setwd(paste(home,'/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/',m.monkey,date,sep=""))
  write.csv(pose_log,file=paste('GROOMLOG_pose_restructured.csv',sep=""),row.names = F)
  write.csv(bodypart_log,file=paste('GROOMLOG_bodypart_restructured.csv',sep=""),row.names = F)

