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
library(tidyverse)

#Load data:
mac =1
if (mac == 1) { home = '~'} else
{home = 'C:/Users/GENERAL'}

setwd(paste(home,'/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output', sep=""))
file = file.choose() # chose the formatted behavior file (xls file)
if (length(grep('partner',file)>0)) {
  if (length(grep('Amos',file))>0) {monkey = "Lovelace"; m.monkey="Amos" } else {monkey = "Sally"; m.monkey="Hooke" }
  date = as.character(substr(file, nchar(file)-33, nchar(file)-23))
} else if (length(grep('neighbor',file)>0)){
  if (length(grep('Amos',file))>0) {monkey = "Amos" } else {monkey = "Hooke"}
  date = as.character(substr(file, nchar(file)-34, nchar(file)-24))
} else { if (length(grep('Amos',file))>0) {monkey = "Amos" } else {monkey = "Hooke"}
  date = as.character(substr(file, nchar(file)-25, nchar(file)-15))
}

log <- xlsx::read.xlsx(file, 1)

#Get all unique behaviors
behavs = as.character(unique(log$Behavior))
behavs = behavs[!(is.na(behavs) | behavs=="")] #remove blank behaviors

#Only keep behaviors, remove session events
to_remove = cbind("Started recording","Stopped recording")
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
  data$start.time=as.numeric(data$start.time); data$end.time=as.numeric(data$end.time)
  data[,c("start.time","end.time")]=data[,c("start.time","end.time")]/1000
  
  new_log=rbind(new_log,data)
}
#behavs[b]

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

#Get block limits
block_order = new_log$Behavior[grep('block',new_log$Behavior)]
block_start = new_log$start.time[grep('block',new_log$Behavior)]
block_end = new_log$end.time[grep('block',new_log$Behavior)]

#Order behaviors
new_log_final = new_log; unique(new_log$Behavior)
new_log_final$Behavior=factor(new_log_final$Behavior,
                              levels=c("Aggression","Proximity", "HIP","Foraging", "Vocalization","SS","Groom Give", "Masturbating","Mounting",
                                       "Submission", "Approach","Yawning","Self-groom","HIS","Other monkeys vocalize", "Lip smack",
                                       "Groom Receive","Leave","Drinking","SP","Pacing/Travel","Scratch","RR", "Butt sniff","Grm prsnt",
                                       "Swinging", "Head Bobbing", "Object Manipulation","watch.neighbor","Camera Sync"))

#Remove NAs (for behavior categories we do not consider here)
new_log_final = new_log_final[!is.na(new_log_final$Behavior),]
max(new_log_final$duration.s)

if (length(which(new_log_final$duration.s<=0))>0){stop("NEGATIVE OR 0 DURATION")}

new_log_final_plot <- new_log_final
new_log_final_plot$Behavior = as.character(new_log_final_plot$Behavior)
new_log_final_plot$Behavior[new_log_final_plot$Behavior=="SS"]="HIS"
new_log_final_plot$Behavior[new_log_final_plot$Behavior=="SP"]="HIP"
new_log_final_plot$Behavior[new_log_final_plot$Behavior=="HIS"]="Threat to subject"
new_log_final_plot$Behavior[new_log_final_plot$Behavior=="HIP"]="Threat to partner"
new_log_final_plot$Behavior[new_log_final_plot$Behavior=="Pacing/Travel"]="Travel"
new_log_final_plot$Behavior[new_log_final_plot$Behavior=="Groom Give"]="Groom partner"
new_log_final_plot$Behavior[new_log_final_plot$Behavior=="Groom Receive"]="Getting groomed"
new_log_final_plot$Behavior[new_log_final_plot$Behavior=="RR"]="Alert"
new_log_final_plot$Behavior[new_log_final_plot$Behavior=="Grm prsnt"]="Limb presentation"
new_log_final_plot$Behavior = factor(new_log_final_plot$Behavior)

new_log_final_plot$Behavior=factor(new_log_final_plot$Behavior,
                              levels=c("Aggression","Proximity", "Threat to partner","Foraging", "Vocalization","Groom partner", "Mounting",
                                      "Approach","Yawning","Self-groom","Threat to subject","Other monkeys vocalize", 
                                       "Getting groomed","Leave","Drinking","Travel","Scratch","Alert", "Limb presentation","Camera Sync"))

#Plot
behavior.log<-ggplot(new_log_final_plot, aes(xmin=start.time, xmax= end.time, ymin=group.min, ymax=group.max))+
  geom_rect(aes(fill=Behavior))+#, colour = "grey50")+
  geom_vline(xintercept = block_end[1])+
  geom_vline(xintercept = block_end[2])+
  #scale_fill_viridis(option="turbo", discrete = TRUE)+
  theme_classic(base_size = 16)+ ylim(0,1)+xlim(0,max(new_log$end.time))+
  xlab('Time since start of recording (in s)')+
  theme(axis.text.y= element_blank(),
        axis.ticks.y = element_blank())#+
#scale_x_continuous(breaks=c(0,600,2000,4000,6000))

#Save plot
if (length(grep('partner',file))>0) {
  setwd(paste(home,'/Dropbox (Penn)/Datalogger/Results/',m.monkey,date,'/Behavior_results/',sep=""))
  ggsave(behavior.log,filename = paste("behavior_log_plot_",monkey,date,".pdf", sep="")) }else
  {setwd(paste(home,'/Dropbox (Penn)/Datalogger/Results/',monkey,date,'/Behavior_results/',sep=""))
    ggsave(behavior.log,filename = paste("behavior_log_plot_",monkey,date,".pdf", sep=""))}

#Add block limits
blocklim = data.frame(matrix(NA, nrow = 3, ncol = ncol(new_log_final))); names(blocklim)=names(new_log_final)
blocklim$Behavior = block_order
blocklim$start.time = block_start
blocklim$end.time = block_end
new_log_final = rbind(new_log_final, blocklim)

#Save to .csv
#output_file = utils::choose.dir(default = "", caption = "Select folder") # choose output directory
# dir <- dirname(output_file)
if (length(grep('partner',file))>0) {
  setwd(paste(home,'/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/',m.monkey,date,sep=""))
  write.csv(new_log_final[,-c(4,5)],file=paste('EVENTLOG_restructured_partner.csv',sep=""),row.names = F)
  } else if (length(grep('neighbor',file))>0) {
    setwd(paste(home,'/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/',monkey,date,sep=""))
    write.csv(new_log_final[,-c(4,5)],file=paste('EVENTLOG_restructured_neighbor.csv',sep=""),row.names = F)
  } else {setwd(paste(home,'/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/',monkey,date,sep=""))
    write.csv(new_log_final[,-c(4,5)],file=paste('EVENTLOG_restructured.csv',sep=""),row.names = F)}

