#Behavioral timeline during pilot recordings of Amos

#Load libraries: 
library(timeline)
library(ggplot2)

#Load data:
setwd('C:/Users/Camille Testard/Desktop/Grants_Fellowships/R37_Project/pilot_data/session3')
behavioral_log = read.csv('behavioral_log_session3.csv')
behavioral_log$Start.Stop=as.character(behavioral_log$Start.Stop)
names(behavioral_log)[1]="Time"

# #For now exclude ambiguous start/stop
# log=behavioral_log[-which(str_length(as.character(behavioral_log$Start.Stop))>5),]
log=behavioral_log

#Reorganize data for plotting
behavs = as.character(unique(log$Behavior))
awakeness = behavs[c(3,4)] #sleep/wake states. Seb forgot awake!
behavs = behavs[-c(1:4,9,10,15)] #behaviors
wakingup_time = log$Time[log$Behavior=="Waking up"];

#For all behaviors outside of sleep/wake cycle:
new_log=data.frame(); b=1
for (b in 1:length(behavs)){
  
  data=log[log$Behavior==behavs[b],c("Time", "Behavior", "Start.Stop")]
  if (length(which(is.na(data$Behavior)))!=0){data=data[-which(is.na(data$Behavior)),]}
  
  interim=split(data,data$Start.Stop)
  names(interim[[1]])[1]="start.time";names(interim[[2]])[1]="stop.time"
  
  data = cbind(interim[[1]],interim[[2]]); data=data[,-c(3,5,6)]
  data=data[,c("Behavior","start.time","stop.time")]
  data[,c("start.time","stop.time")]=data[,c("start.time","stop.time")]
  
  new_log=rbind(new_log,data)
}

#Add sleeping:
times = c(log$Time[which(log$Behavior=="Waking up")])
new_log[nrow(new_log)+1,"Behavior"]="Asleep"
new_log[nrow(new_log),"start.time"]=1; new_log[nrow(new_log),"stop.time"]=times[1];
new_log$duration.s = new_log$stop.time-new_log$start.time

if (length(which(new_log$duration.s<0)) !=0){print('ERROR NEGATIVE DURATION')}

new_log$group.min=0; new_log$group.max=1
new_log$Behavior=as.character(new_log$Behavior)

#Add resting:
new_log = new_log[order(new_log$start.time),]
start.times=new_log$start.time[2:length(new_log$start.time)]
stop.times=new_log$stop.time[1:length(new_log$stop.time)-1]
intervals = as.numeric(start.times-stop.times)
free_intervals = which(intervals>0)

resting=data.frame()
for (i in 1:length(free_intervals)){
  resting[i,'start.time']=stop.times[free_intervals[i]]
  resting[i,'stop.time']=start.times[free_intervals[i]]
}
resting$duration.s = resting$stop.time-resting$start.time
resting$Behavior="Resting";
resting$group.min=0; resting$group.max=1; resting=resting[,names(new_log)]

new_log=rbind(new_log,resting)
new_log = new_log[order(new_log$start.time),]
new_log$Behavior=factor(new_log$Behavior,
                        levels=c("Aggression","Proximity","Grooming (giving)", "Grooming (receiving)",
                                 "Feeding", "Submission","Yawning",
                                 "Resting","Walking","Drinking","Asleep","Self-directed behavior"))
# #Add wakefulness state
# new_log$wake.state = 'awake'
# new_log$wake.state[1] = 'sleeping'
# wakingup_idx = max(which(new_log$start.time<fully_awake_time))
# new_log$wake.state[2:wakingup_idx] = 'wakingup'
# stress_idx = which(new_log$start.time>=stress_time)
# new_log$wake.state[stress_idx] = 'stressed'

#Save to .csv
write.csv(new_log[,-c(5,6)],file='behav_log_restructured.csv',row.names = F)

