#Get_freq_prop_behav.R

#Load libraries: 
library(ggplot2) 
library(utils)
library(xlsx)

#Load data:
sessions = c("Amos_2021-07-29","Hooke_2021-08-02","Amos_2021-08-03","Hooke_2021-08-12")
session_length = c(8897, 7828, 8897, 7184) #Check in OFS files!

all_logs = data.frame()
for (s in 1:length(sessions)){
  setwd(paste('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/', sessions[s],sep=""))
  log = read.csv(paste('EVENTLOG_restructured_',sessions[s],'.csv',sep=""))
  log$session_num = s
  all_logs= rbind(all_logs, log)
}

#Rename certain behaviors:
all_logs$Behavior[all_logs$Behavior=="RR"]= "Rowdy room"
all_logs$Behavior[all_logs$Behavior=="Pacing/Travel"]= "Travel"
all_logs$Behavior[all_logs$Behavior=="HIS"]= "Human Intrusion Subject"
all_logs$Behavior[all_logs$Behavior=="HIP"]= "Human Intrusion Partner"
all_logs$Behavior[all_logs$Behavior=="SS"]= "Squeeze Subject"
all_logs$Behavior[all_logs$Behavior=="SP"]= "Squeeze Partner"

#Get frequency table per behavior
behavs = unique(all_logs$Behavior);
behav_categ = c("social","social","asocial","social","asocial", "social",
                "social","social","asocial","asocial","social","asocial",
                "asocial","social","social","asocial","social","social","social","social","asocial")
cbind(behavs, behav_categ)
behav_freq_table = data.frame(); b=1; s=1
behav_freq_table_full = data.frame()
for (s in 1:length(sessions)){
  for (b in 1:length(behavs)){
    behav_freq_table[b,'Behavior'] = behavs[b]
    behav_freq_table[b,'categ'] = behav_categ[b]
    behav_freq_table[b,'session'] = s
    behav_freq_table[b,'obs.time'] = session_length[s]
    idx = which(all_logs$Behavior == behavs[b] & all_logs$session_num == s)
    behav_freq_table[b,'duration'] = sum(all_logs$duration.s[idx])
    behav_freq_table[b,'proportion'] = behav_freq_table[b,'duration']/behav_freq_table[b,'obs.time']
    behav_freq_table[b,'num.events'] = length(idx)
  }
 
  setwd(paste('~/Dropbox (Penn)/Datalogger/Results/', sessions[s],"/",sep=""))
  
  p.dur<- ggplot(data=behav_freq_table, aes(x=reorder(Behavior, duration), y=duration, fill=categ)) +
    geom_bar(stat="identity")+ xlab('Behaviors')+ ylab('Duration (in s)')+ 
    labs(fill='Category')+ 
    coord_flip()+ theme_classic(base_size = 18) + ggtitle(sessions[s])
  ggsave(paste("Behavior_duration_", sessions[s],".png",sep=""))
  
  p.prop<- ggplot(data=behav_freq_table, aes(x=reorder(Behavior, proportion), y=proportion, fill=categ)) +
    geom_bar(stat="identity")+ xlab('Behaviors')+ ylab('Proportion')+ 
    labs(fill='Category')+ 
    coord_flip()+ theme_classic(base_size = 18) + ggtitle(sessions[s])
  ggsave(paste("Behavior_proportion_", sessions[s],".png",sep=""))
  
  p.frq<- ggplot(data=behav_freq_table, aes(x=reorder(Behavior, num.events), y=num.events, fill=categ)) +
    geom_bar(stat="identity")+ xlab('Behaviors')+ ylab('# Occurences')+
    coord_flip()+ theme_classic(base_size = 18)+ ggtitle(sessions[s])
  ggsave(paste("Behavior_frequency_", sessions[s],".png",sep=""))
  
  behav_freq_table_full = rbind(behav_freq_table_full, behav_freq_table)
}

setwd('~/Dropbox (Penn)/Datalogger/Results/')
ggplot(data=behav_freq_table_full, aes(x=reorder(Behavior, duration), y=duration, fill=categ)) +
  geom_bar(stat="identity")+ xlab('Behaviors')+ ylab('Duration (in s)')+ 
  labs(fill='Category')+ 
  coord_flip()+ theme_classic(base_size = 18) + ggtitle("All session cumulated")
ggsave("Behavior_proportion_allSessions.png")


ggplot(data=behav_freq_table_full, aes(x=reorder(Behavior, num.events), y=num.events, fill=categ)) +
  geom_bar(stat="identity")+ xlab('Behaviors')+ ylab('# Occurences')+
  coord_flip()+ theme_classic(base_size = 18)+ ggtitle("All session cumulated")
ggsave(paste("Behavior_frequency_", sessions[s],".png",sep=""))

#Get proportion 