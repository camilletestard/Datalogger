#Format_behavior_log.R
#This script re-formats the deuteron behavior log in a way that is easily interpretable for 
#neural data analysis. It also plots behavior during a session. 
# Certain red flags to keep in mind: 
#Camille Testard - November 2021

#Load libraries: 
library(timeline)
library(ggplot2) 
library(utils)
library(xlsx)

#Load data:
file = file.choose() # chose the formatted behavior file
monkey = "Amos"
log = read.xlsx(file, sheetIndex = 1)
#log = read.csv(file)

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
for (b in 1:length(behavs)){ #For all behaviors
  
  data=log[log$Behavior==behavs[b],c("Time", "Behavior", "Start.end")] #Only keep useful columns
  if (length(which(is.na(data$Behavior)))!=0){data=data[-which(is.na(data$Behavior)),]}# remove rows with empty/NA behavior
  
  interim=split(data,data$Start.end) #Split the data by behavior start and behavior end
  names(interim[[2]])[1]="start.time";names(interim[[1]])[1]="end.time" #rename columns
  
  data = cbind(interim[[1]],interim[[2]]); data=data[,-c(3,5,6)]
  data=data[,c("Behavior","start.time","end.time")]
  data[,c("start.time","end.time")]=data[,c("start.time","end.time")]/1000
  
  new_log=rbind(new_log,data)
}

#Note: The bit of code above will have an error if there aren't the same
#number of start and end for any one behavior. When an error occurs, check for which behavior it occurs
# and adjust the excel file accordingly (usually an end is missing, or there are 2 starts, use the videos to check). 
# The script stops at the behavior that has an issue, number "b". You can check in the character vector
# 'behavs' which behavior causes issues.

#Format for plotting:
new_log$group.min=0; new_log$group.max=1
new_log$Behavior=as.character(new_log$Behavior)

#Order new log chronologically
new_log = new_log[order(new_log$start.time),]
new_log$duration.s = new_log$end.time - new_log$start.time

#Find block limits
block_limits = new_log$start.time[c(which(new_log$Behavior=="Unpair"), grep('neighbor',new_log$Behavior))]
block_limits = block_limits[order(block_limits)];

#Order behaviors
new_log_final = new_log; unique(new_log$Behavior)
new_log_final$Behavior=factor(new_log_final$Behavior, 
                              levels=c("Aggression","Proximity","Groom Give", "HIP","Foraging", "Vocalization","SS", "Masturbating",
                                       "Submission", "Approach","Yawning","Self-groom","HIS","Other monkeys vocalize",
                                       "Groom Receive","Leave","Drinking","SP","Pacing/Travel","Scratch","RR"))

#Remove NAs (for behavior categories we do not consider here)
new_log_final = new_log_final[!is.na(new_log_final$Behavior),]

if (length(which(new_log_final$duration.s<0))>0){stop("NEGATIVE DURATION")}

#Plot
behavior.log<-ggplot(new_log_final, aes(xmin=start.time, xmax= end.time, ymin=group.min, ymax=group.max))+
  geom_rect(aes(fill=Behavior))+#, colour = "grey50")+
  theme_classic(base_size = 20)+ ylim(0,1)+xlim(0,max(new_log$end.time))+
  xlab('Time since start of recording (in s)')+
  theme(axis.text.y= element_blank(),
        axis.ticks.y = element_blank())#+
#scale_x_continuous(breaks=c(0,600,2000,4000,6000))

#Save plot
ggsave(behavior.log,filename = paste("behavior_log_plot_",monkey,as.character(substr(file, 103, 113)),".png", sep=""))

#Add block limits
blocklim = data.frame(matrix(NA, nrow = 3, ncol = ncol(new_log_final))); names(blocklim)=names(new_log_final)
blocklim$Behavior[1]='Pair1'; 
blocklim$start.time[1]=1; 
blocklim$end.time[1]=round(block_limits[1]); 

blocklim$Behavior[2]='Pair2'; 
blocklim$start.time[2]=round(block_limits[1]); 
blocklim$end.time[2]=round(block_limits[2]); 

blocklim$Behavior[3]='Alone'; 
blocklim$start.time[3]=round(block_limits[2]); 
blocklim$end.time[3]=round(log$Time[nrow(log)]/1000); 

new_log_final = rbind(new_log_final, blocklim)

#Save to .csv
#output_file = utils::choose.dir(default = "", caption = "Select folder") # choose output directory
# dir <- dirname(output_file)
setwd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
write.csv(new_log_final[,-c(4,5)],file=paste('EVENTLOG_restructured_',monkey,as.character(substr(file, 98, 108)),'.csv',sep=""),row.names = F)

