#Get frequency of behaviors for pilot data

duration_session1 = 7606;
duration_session3 = 7466;
duration_session4 = 9119;
total.dur.s = duration_session1+duration_session3+duration_session4
total.dur.h = total.dur.s/3600

#Load data:
setwd('C:/Users/Camille Testard/Desktop/Grants_Fellowships/R37_Project/pilot_data/all_behavior_sessions')
behavioral_log1 = read.csv('behav_log_restructured_session1.csv')
behavioral_log3 = read.csv('behav_log_restructured_session3.csv')
behavioral_log4 = read.csv('behav_log_restructured_session4.csv')

#Get columns of interest
behavioral_log1 = behavioral_log1[,c("Behavior","duration.s")]
behavioral_log3 = behavioral_log3[,c("Behavior","duration.s")]
behavioral_log4 = behavioral_log4[,c("Behavior","duration.s")]

#Combine data
all_logs = rbind(behavioral_log1,behavioral_log3,behavioral_log4)

#Combine grooming to facilitate
idx = all_logs$Behavior=="Grooming (receiving)" | all_logs$Behavior=="Grooming (giving)" 
all_logs$Behavior[idx] = "Grooming"
all_logs$Behavior[all_logs$Behavior=="Asleep"] = "Sleeping"

#Get frequency table per behavior
behavs = unique(all_logs$Behavior);
behav_freq_table = data.frame(); b=1
for (b in 1:length(behavs)){
  behav_freq_table[b,'Behavior'] = behavs[b]
  idx = which(all_logs$Behavior == behavs[b])
  behav_freq_table[b,'duration/hrs observed'] = sum(all_logs$duration.s[idx])/total.dur.h
  behav_freq_table[b,'#occurrence/hrs observed'] = length(idx)/total.dur.h
}

#Save to .csv
write.csv( behav_freq_table,file=' behav_freq_table.csv',row.names = F)