#Get_freq_prop_behav.R

#Load libraries: 
library(ggplot2) 
library(utils)
library(xlsx)
library(dplyr)
library(plotrix)

#Load data:
sessions = c("Amos_2021-07-29","Hooke_2021-08-02","Amos_2021-08-03","Hooke_2021-08-05",
             "Amos_2021-08-09","Hooke_2021-08-12","Amos_2021-08-13","Amos_2021-08-16",
             "Hooke_2021-08-19","Amos_2021-08-20","Hooke_2021-08-21","Hooke_2021-08-29",
             "Hooke_2021-09-03")
a_sessions = c(1,3,5,7,8,10)
h_sessions = c(2,4,6,9,11,12,13)

all_logs = data.frame(); session_length = vector(); s=1
for (s in 1:length(sessions)){
  setwd(paste('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/', sessions[s],sep=""))
  log = read.csv(paste('EVENTLOG_restructured.csv',sep=""))
  log$session_num = s
  session_length[s] = log$end.time[nrow(log)]; #This should be correct within a few seconds.
  all_logs= rbind(all_logs, log)
}

#Rename certain behaviors:
all_logs$Behavior[all_logs$Behavior=="RR"]= "Rowdy room"
all_logs$Behavior[all_logs$Behavior=="Groom Give"]= "Groom partner"
all_logs$Behavior[all_logs$Behavior=="Groom Receive"]= "Getting groomed"
all_logs$Behavior[all_logs$Behavior=="Pacing/Travel"]= "Travel"
all_logs$Behavior[all_logs$Behavior=="HIS"]= "Threat to subject"
all_logs$Behavior[all_logs$Behavior=="HIP"]= "Threat to partner"
all_logs$Behavior[all_logs$Behavior=="SS"]= "Threat to subject"
all_logs$Behavior[all_logs$Behavior=="SP"]= "Threat to partner"
all_logs$Behavior[all_logs$Behavior=="Grm prsnt"]= "Groom solicitation"

#Remove blocks
all_logs = all_logs[-grep("block",all_logs$Behavior),]

#Get frequency table per behavior
behavs = unique(all_logs$Behavior);
behav_categ_social = c("social","social","asocial","social","asocial", "social",
                "social","social","social","social","asocial","asocial",
                "social","social","asocial","asocial","social","social","asocial",
                "social","social","asocial", "social")
behav_categ_valence = c("affiliative", "affiliative","neutral","neutral","neutral",
                        "agonistic","neutral","affiliative","affiliative", "affiliative",
                        "neutral","neutral","affiliative","neutral","neutral","neutral",
                        "agonistic","agonistic","neutral","affiliative",
                        "agonistic","neutral","neutral")
cbind(behavs, behav_categ_social,behav_categ_valence)
behav_freq_table = data.frame(); b=1; s=1
behav_freq_table_full = data.frame()

for (s in 1:length(sessions)){
  
  for (b in 1:length(behavs)){
    behav_freq_table[b,'Behavior'] = behavs[b]
    behav_freq_table[b,'social_categ'] = behav_categ_social[b]
    behav_freq_table[b,'valence_categ'] = behav_categ_valence[b]
    behav_freq_table[b,'session'] = s
    behav_freq_table[b,'obs.time'] = session_length[s]
    behav_freq_table[b,'total.events'] = nrow(all_logs)
    idx = which(all_logs$Behavior == behavs[b] & all_logs$session_num == s)
    behav_freq_table[b,'duration'] = sum(all_logs$duration.s[idx])
    behav_freq_table[b,'proportion'] = behav_freq_table[b,'duration']/behav_freq_table[b,'obs.time']
    behav_freq_table[b,'num.events'] = length(idx)
  }
 
  setwd(paste('~/Dropbox (Penn)/Datalogger/Results/', sessions[s],"/Behavior_results",sep=""))
  
  #Pie chart based on duration
  png(file=paste("Behavior_duration_", sessions[s],".png",sep=""))
  data = behav_freq_table[intersect(which(behav_freq_table$duration>0), which(behav_freq_table$Behavior!= "Proximity")),]
  slices <- data$duration
  lbls <- data$Behavior
  pie(slices, labels = lbls, main="Behaviors observed (duration)", cex=0.6)
 # pie3D(slices, labels = lbls, explode=0.1, main="Behaviors observed", labelcex=0.3)
  dev.off()
  
  #Pie chart based on number of events:
  png(file=paste("Behavior_events_", sessions[s],".png",sep=""))
  slices <- data$num.events
  lbls <- data$Behavior
  pie(slices, labels = lbls, main="Behaviors observed (events)", cex=0.6)
  dev.off()
  
p.dur<- ggplot(data=behav_freq_table, aes(x=reorder(Behavior, duration), y=duration, fill=categ)) +
    geom_bar(stat="identity")+ coord_polar("y", start=0)+
    xlab('Behaviors')+ ylab('Duration (in s)')+ 
    labs(fill='Category')+ 
    coord_flip()+ theme_classic(base_size = 14) + ggtitle(sessions[s])
  #ggsave(paste("Behavior_duration_", sessions[s],".png",sep=""))
  
  p.frq<- ggplot(data=behav_freq_table, aes(x=reorder(Behavior, num.events), y=num.events, fill=categ)) +
    geom_bar(stat="identity")+ xlab('Behaviors')+ ylab('# Occurences')+
    coord_flip()+ theme_classic(base_size = 18)+ ggtitle(sessions[s])
  #ggsave(paste("Behavior_frequency_", sessions[s],".png",sep=""))
  
  behav_freq_table_full = rbind(behav_freq_table_full, behav_freq_table)
}

setwd('~/Dropbox (Penn)/Datalogger/Results/All_sessions/')


for (b in 1:length(behavs)){
  behav_freq_table[b,'Behavior'] = behavs[b]
  behav_freq_table[b,'social_categ'] = behav_categ_social[b]
  behav_freq_table[b,'valence_categ'] = behav_categ_valence[b]
  behav_freq_table[b,'session'] = s
  behav_freq_table[b,'obs.time'] = session_length[s]
  behav_freq_table[b,'total.events'] = nrow(all_logs)
  idx = which(all_logs$Behavior == behavs[b])
  behav_freq_table[b,'duration'] = sum(all_logs$duration.s[idx])
  behav_freq_table[b,'proportion'] = behav_freq_table[b,'duration']/behav_freq_table[b,'obs.time']
  behav_freq_table[b,'num.events'] = length(idx)
}

#eps(file=paste("Behavior_duration_allSessions.png",sep=""))
data = behav_freq_table[intersect(which(behav_freq_table$duration>100), which(behav_freq_table$Behavior!= "Proximity")),]
slices <- data$duration
lbls <- data$Behavior
pie(slices, labels = lbls, main="Behaviors observed (duration)", cex=1, col=rainbow(nrow(data)))
# pie3D(slices, labels = lbls, explode=0.1, main="Behaviors observed", labelcex=0.3)
#dev.off()

data2 = data.frame()
data2[1,"social_categ"] = "social"
data2[1,"duration"] = sum(data$duration[data$social_categ=="social"]) 
data2[2,"social_categ"] = "asocial"
data2[2,"duration"] = sum(data$duration[data$social_categ=="asocial"]) 
slices <- data2$duration
lbls <- data2$social_categ
pie(slices, labels = lbls, main="Behaviors observed (duration)", cex=2, col=cm.colors(2))
data2$duration[1]/sum(data2$duration)

data3 = data.frame()
data3[1,"valence"] = "affiliative"
data3[1,"duration"] = sum(data$duration[data$valence_categ=="affiliative"]) 
data3[2,"valence"] = "agonistic"
data3[2,"duration"] = sum(data$duration[data$valence_categ=="agonistic"]) 
slices <- data3$duration
lbls <- data3$valence
pie(slices, labels = lbls, main="Behaviors observed (duration)", cex=2, col=heat.colors(2))
data3$duration[1]/sum(data3$duration)

#Format behavior freq table
behav_freq_table_full<-behav_freq_table_full[behav_freq_table_full$Behavior!="Proximity",]
behav_freq_table_full$duration[is.na(behav_freq_table_full$duration)]=0

ggplot(data=behav_freq_table_full, aes(x=reorder(Behavior, num.events), y=num.events, fill=social_categ)) +
  geom_bar(stat="identity")+ xlab('Behaviors')+ ylab('# Occurences')+
  coord_flip()+ theme_classic(base_size = 18)+ ggtitle("All session cumulated")
ggsave(paste("Behavior_frequency_allSessions.png",sep=""))

ggplot(data=behav_freq_table_full, aes(x=reorder(Behavior, duration), y=duration)) +
  geom_bar(stat="identity")+ xlab('Behaviors')+ ylab('Duration (s)')+
  coord_flip()+ theme_classic(base_size = 16)+ ggtitle("All session cumulated")


ggplot(data=behav_freq_table_full, aes(x=reorder(Behavior, proportion), y=proportion/4, fill=categ)) +
  geom_bar(stat="identity")+ xlab('Behaviors')+ ylab('Proportion')+
  coord_flip()+ theme_classic(base_size = 18)+ ggtitle("All session cumulated")
ggsave(paste("Behavior_proportion_allSessions.png",sep=""))
