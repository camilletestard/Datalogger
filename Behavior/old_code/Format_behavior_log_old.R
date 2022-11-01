#Format_behavior_log.R
#This script re-formats the deuteron behavior log in a way that is easily interpretable for neural data analysis

#Load libraries: 
library(timeline)
library(ggplot2) 
library(utils)

setwd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/EventLogs/EventLogs_Done')

#Load data:
file = file.choose() # chose the formatted behavior file
behavioral_log = read.csv(file)

# #For now exclude ambiguous start/end
# log=behavioral_log[-which(str_length(as.character(behavioral_log$Start.end))>5),]
log=behavioral_log

#Reorganize data for plotting
behavs = as.character(unique(log$Behavior))
behavs = behavs[!(is.na(behavs) | behavs=="")]
#awakeness = behavs[c(2,3,5,11)] #sleep/wake states

#Only keep behaviors. STOP!!!! need to check for each session because behavior ordering differs from one session to the next.
#behavs = behavs[-c(1,2)] #only keep behaviors
to_remove = cbind("Started recording", "Camera Sync", "Stopped recording")
if (any(is.na(match(to_remove, behavs)))) {
  behavs
  to_remove = to_remove[!is.na(match(to_remove, behavs))]
}
behavs = behavs[-match(to_remove, behavs)]  #behavs[-c(1,2,length(behavs))] #only keep behaviors

#For all behaviors outside of sleep/wake cycle:
new_log=data.frame(); b=1
for (b in 11:length(behavs)){
  
  data=log[log$Behavior==behavs[b],c("Time", "Behavior", "Start.end")]
  if (length(which(is.na(data$Behavior)))!=0){data=data[-which(is.na(data$Behavior)),]}# remove rows with empty/NA behavior
  
  interim=split(data,data$Start.end)
  names(interim[[2]])[1]="start.time";names(interim[[1]])[1]="end.time"
  
  data = cbind(interim[[1]],interim[[2]]); data=data[,-c(3,5,6)]
  data=data[,c("Behavior","start.time","end.time")]
  data[,c("start.time","end.time")]=data[,c("start.time","end.time")]/1000
  
  new_log=rbind(new_log,data)
}

# #Add sleeping:
# times = c(log$Time[which(log$Behavior=="wakingup")]/1000,
#           log$Time[which(log$Behavior=="awake")]/1000)
# new_log[nrow(new_log)+1,"Behavior"]="Sleeping"
# new_log[nrow(new_log),"start.time"]=1; new_log[nrow(new_log),"end.time"]=times[1];
# new_log$duration.s = new_log$end.time-new_log$start.time

new_log$group.min=0; new_log$group.max=1
new_log$Behavior=as.character(new_log$Behavior)

# #Add resting (i.e. whenever no behavior was recorded):
# new_log = new_log[order(new_log$start.time),]
# start.times=new_log$start.time[2:length(new_log$start.time)]
# end.times=new_log$end.time[1:length(new_log$end.time)-1]
# intervals = as.numeric(start.times-end.times)
# free_intervals = which(intervals>0)
# 
# resting=data.frame()
# for (i in 1:length(free_intervals)){
#   resting[i,'start.time']=end.times[free_intervals[i]]
#   resting[i,'end.time']=start.times[free_intervals[i]]
# }
# resting$duration.s = resting$end.time-resting$start.time
# resting$Behavior="Resting"; resting$group.min=0; resting$group.max=1; resting=resting[,names(new_log)]
# 
# new_log=rbind(new_log,resting)

#Order new log chronologically
new_log = new_log[order(new_log$start.time),]
new_log$duration.s = new_log$end.time - new_log$start.time

if (length(which(new_log$duration.s<0))>0){stop("NEGATIVE DURATION")}

#Order behaviors for later plotting
new_log_final = new_log; unique(new_log$Behavior)
new_log_final$Behavior=factor(new_log_final$Behavior, 
                        levels=c("Aggression","Proximity","Groom Give", "HIP","Foraging", "Vocalization","SS", "Masturbating",
                                 "Submission", "Approach","Yawning","Self-groom","HIS","Other monkeys vocalize",
                                 "Groom Receive","Leave","Drinking","SP","Pacing/Travel","Scratch","RR"))

#Remove NAs (for behavior categories we do not consider here)
new_log_final = new_log_final[!is.na(new_log_final$Behavior),]

# levels=c("Aggression","Proximity","Grooming", 
#          "Feeding", "Submission",
#          "Resting","Walking","Drinking","Sleeping","Self-directed behavior"))

# #Add wakefulness state
# new_log$wake.state = 'awake'
# new_log$wake.state[1] = 'sleeping'
# new_log$wake.state[2:11] = 'wakingup'

#Save to .csv
# output_file = utils::choose.dir(default = "", caption = "Select folder") # choose output directory
# dir <- dirname(output_file)
setwd('~/Dropbox (Penn)/Deuteron_Backup/Deuteron_Data_Backup/Ready to analyze output/')
write.csv(new_log_final[,-c(4,5)],file=paste('EVENTLOG_restructured',as.character(substr(file, 95, 113)),'.csv',sep=""),row.names = F)

#Plot
behavior.log<-ggplot(new_log_final, aes(xmin=start.time, xmax= end.time, ymin=group.min, ymax=group.max))+
  geom_rect(aes(fill=Behavior))+#, colour = "grey50")+
  theme_classic(base_size = 20)+ ylim(0,1)+xlim(0,max(new_log$end.time))+
  xlab('Time since start of recording (in s)')+
  theme(axis.text.y= element_blank(),
        axis.ticks.y = element_blank())#+
  #scale_x_continuous(breaks=c(0,600,2000,4000,6000))

ggsave(behavior.log,filename = paste("behavior_log_plot",as.character(substr(file, 95, 113)),".png"))

# #Add sleep/wake sates:
# new_log2=data.frame()
# times = c(log$Time[which(log$Behavior=="sleep")]/1000,
#           log$Time[which(log$Behavior=="wakingup")]/1000,
#           log$Time[which(log$Behavior=="awake")]/1000)
# new_log2[1,c("Behavior","start.time","end.time")]=c("sleep",1,times[2])
# new_log2[2,c("Behavior","start.time","end.time")]=c("waking.up",times[2],times[3])
# new_log2[3,c("Behavior","start.time","end.time")]=c("awake",times[3],max(new_log$end.time))
# new_log2$start.time = as.numeric(new_log2$start.time); new_log2$end.time = as.numeric(new_log2$end.time)
# new_log2$duration.s = new_log2$end.time-new_log2$start.time
# new_log2$group.min=1; new_log2$group.max=2

# #Combine the two
# final_log=rbind(new_log,new_log2)

# ggplot(final_log, aes(xmin=start.time, xmax= end.time, ymin=group.min, ymax=group.max))+
#   geom_rect(aes(fill=Behavior))+
#   theme_classic(base_size = 20)+
#   xlab('Time since start of recording (in s)')+
#   scale_y_discrete(limits=c(0.5, 1.5), labels=c("Behavior","Sleep/Wake State"))+
#   theme(axis.text.y= element_text(angle=45))

#other package:
# timeline(new_log, label.col='Behavior', start.col='start.time', end.col='end.time',
#          group.col='Group')
