#Format_behavior_log.R
#This script re-formats the deuteron behavior log in a way that is easily interpretable for
#neural data analysis. It also plots behavior during a session.
# Certain red flags to keep in mind:
#Camille Testard - September 2021
                  
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

# #Split all behavior events by start and end:
# new_log=data.frame(); b=1
# for (b in 1:length(behavs)){ #For all behaviors
#   
#     data=log[log$Behavior==behavs[b],c("Time", "Behavior", "Start.end")] #Only keep useful columns
#   if (length(which(is.na(data$Behavior)))!=0){data=data[-which(is.na(data$Behavior)),]}# remove rows with empty/NA behavior
#   
#   interim=split(data,data$Start.end) #Split the data by behavior start and behavior end
#   names(interim[[2]])[1]="start.time";names(interim[[1]])[1]="end.time" #rename columns
#   
#   data = cbind(interim[[1]],interim[[2]]); data=data[,-c(3,5,6)]
#   data=data[,c("Behavior","start.time","end.time")]
#   data$start.time=as.numeric(data$start.time); data$end.time=as.numeric(data$end.time)
#   data[,c("start.time","end.time")]=data[,c("start.time","end.time")]/1000
#   
#   new_log=rbind(new_log,data)
# }
# #behavs[b]
# 
# #Note: The bit of code above will have an error if there aren't the same
# #number of start and end for any one behavior. When an error occurs, check for which behavior it occurs
# # and adjust the excel file accordingly (usually an end is missing, or there are 2 starts, use the videos to check).
# # The script stops at the behavior that has an issue, number "b". You can check in the character vector
# # 'behavs' which behavior causes issues.
# 
# #Format for plotting:
# new_log$group.min=0; new_log$group.max=1
# new_log$Behavior=as.character(new_log$Behavior)
# 
# #Order new log chronologically
# new_log = new_log[order(new_log$start.time),]
# new_log$duration.s = new_log$end.time - new_log$start.time
# 
# #Get block limits
# block_order = new_log$Behavior[grep('block',new_log$Behavior)]
# block_start = new_log$start.time[grep('block',new_log$Behavior)]
# block_end = new_log$end.time[grep('block',new_log$Behavior)]
# 
# #Order behaviors
# new_log_final = new_log; unique(new_log$Behavior)
# new_log_final$Behavior=factor(new_log_final$Behavior,
#                               levels=c("Groom Give","Groom Receive"))
# 
# if (length(which(new_log_final$duration.s<=0))>0){stop("NEGATIVE OR 0 DURATION")}
# 
# #Add block limits
# blocklim = data.frame(matrix(NA, nrow = 3, ncol = ncol(new_log_final))); names(blocklim)=names(new_log_final)
# blocklim$Behavior = block_order
# blocklim$start.time = block_start
# blocklim$end.time = block_end
# new_log_final = rbind(new_log_final, blocklim)

#Save to .csv
#output_file = utils::choose.dir(default = "", caption = "Select folder") # choose output directory
# dir <- dirname(output_file)
  setwd(paste(home,'/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/',m.monkey,date,sep=""))
  write.csv(pose_log,file=paste('GROOMLOG_pose_restructured.csv',sep=""),row.names = F)
  write.csv(bodypart_log,file=paste('GROOMLOG_bodypart_restructured.csv',sep=""),row.names = F)

