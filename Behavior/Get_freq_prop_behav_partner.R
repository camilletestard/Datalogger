#Get_freq_prop_behav.R

#Load libraries: 
library(ggplot2) 
library(utils)
library(xlsx)
library(dplyr)
library(plotrix)

#Load data:
sessions = c("Amos_2021-07-29","Hooke_2021-08-02","Amos_2021-08-03",
             "Hooke_2021-08-05","Amos_2021-08-09","Hooke_2021-08-12")
a_sessions = c(1,3,5)
h_sessions = c(2,4,6)

setwd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
behav_categs = read.csv(file="behav_categs.csv")

all_logs_subject = data.frame(); all_logs_partner = data.frame(); session_length = vector(); s=4
for (s in 1:length(sessions)){
  setwd(paste('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/', sessions[s],sep=""))
  log_subject = read.csv(paste('EVENTLOG_restructured.csv',sep=""))
  log_subject$session_num = s
  
  log_partner = read.csv(paste('EVENTLOG_restructured_partner.csv',sep=""))
  log_partner$session_num = s
  
  session_length[s] = log_subject$end.time[nrow(log_subject)]; #This should be correct within a few seconds.
  
  all_logs_subject= rbind(all_logs_subject, log_subject)
  all_logs_partner= rbind(all_logs_partner, log_partner)
}

#Rename certain behaviors:
all_logs_subject$Behavior[all_logs_subject$Behavior=="RR"]= "Rowdy room"
all_logs_subject$Behavior[all_logs_subject$Behavior=="Groom Give"]= "Groom partner"
all_logs_subject$Behavior[all_logs_subject$Behavior=="Groom Receive"]= "Getting groomed"
all_logs_subject$Behavior[all_logs_subject$Behavior=="Pacing/Travel"]= "Travel"
all_logs_subject$Behavior[all_logs_subject$Behavior=="HIS"]= "Threat to subject"
all_logs_subject$Behavior[all_logs_subject$Behavior=="HIP"]= "Threat to partner"
all_logs_subject$Behavior[all_logs_subject$Behavior=="SS"]= "Threat to subject"
all_logs_subject$Behavior[all_logs_subject$Behavior=="SP"]= "Threat to partner"
all_logs_subject$Behavior[all_logs_subject$Behavior=="Grm prsnt"]= "Limb presentation"

all_logs_partner$Behavior[all_logs_partner$Behavior=="RR"]= "Rowdy room"
all_logs_partner$Behavior[all_logs_partner$Behavior=="Groom Give"]= "Groom partner"
all_logs_partner$Behavior[all_logs_partner$Behavior=="Groom Receive"]= "Getting groomed"
all_logs_partner$Behavior[all_logs_partner$Behavior=="Pacing/Travel"]= "Travel"
all_logs_partner$Behavior[all_logs_partner$Behavior=="HIS"]= "Threat to subject"
all_logs_partner$Behavior[all_logs_partner$Behavior=="HIP"]= "Threat to partner"
all_logs_partner$Behavior[all_logs_partner$Behavior=="SS"]= "Threat to subject"
all_logs_partner$Behavior[all_logs_partner$Behavior=="SP"]= "Threat to partner"
all_logs_partner$Behavior[all_logs_partner$Behavior=="Grm prsnt"]= "Limb presentation"

#Remove blocks
all_logs_subject = all_logs_subject[-grep("block",all_logs_subject$Behavior),]
all_logs_partner = all_logs_partner[-grep("block",all_logs_partner$Behavior),]

#Remove NAs
which(is.na(all_logs_subject), arr.ind = T)
which(is.na(all_logs_partner), arr.ind = T)
all_logs_subject=all_logs_subject[!is.na(all_logs_subject$Behavior),]
all_logs_partner=all_logs_partner[!is.na(all_logs_partner$Behavior),]

#Get behavior category labels
behavs =behav_categs$behavs
# behavs = unique(c(all_logs_partner$Behavior, all_logs_subject$Behavior));
# behav_categ_social = c("social","social","asocial","asocial","social","asocial", "social",
#                        "social","asocial","social","asocial","asocial","social",
#                        "asocial","asocial","social","social","social","asocial","asocial",
#                        "asocial","social","asocial", "social", "social", "social")
# behav_categ_valence = c("affiliative", "neutral","neutral","neutral","affiliative",
#                         "neutral","agonistic","affiliative","neutral", "affiliative",
#                         "neutral","neutral","affiliative","neutral","neutral","agonistic",
#                         "agonistic","affiliative","neutral",
#                         "neutral","neutral","affiliative","neutral","affiliative",
#                         "agonistic","neutral")
# behav_categ_recip = c("reciprocal","Not", "Not","Not","Not","reciprocal","Not","reciprocal",
#                         "Not","reciprocal","Not","Not", "reciprocal",
#                         "Not","Not","reciprocal","reciprocal","Not","reciprocal","Not","Not",
#                         "Not","Not","Not","Not","Not")
# behav_categs = cbind(behavs, behav_categ_social,behav_categ_valence, behav_categ_recip)
# write.csv(file="behav_categs.csv",behav_categs, row.names = F)


#########################################################
## Get behavior partitioning for all sessions

setwd('~/Dropbox (Penn)/Datalogger/Results/All_sessions/')

behav_freq_table = data.frame(); behav_freq_table_all = data.frame(); b=11;
for (s in 1:length(sessions)){
  for (b in 1:length(behavs)){
    
    behav_freq_table[b,'Behavior'] = behavs[b]
    behav_freq_table[b,'social_categ'] = behav_categs$behav_categ_social[b]
    behav_freq_table[b,'valence_categ'] = behav_categs$behav_categ_valence[b]
    behav_freq_table[b,'recip_categ'] = behav_categs$behav_categ_recip[b]
    behav_freq_table[b,'session']=s
    
    idx_partner = which(all_logs_partner$Behavior == behavs[b] & all_logs_partner$session_num==s)
    sec_partner = vector()
    if (length(idx_partner)>0) {
      for (inter in 1:length(idx_partner)){
        sec_partner = c(sec_partner, round(c(all_logs_partner$start.time[idx_partner[inter]]:
                                               (all_logs_partner$end.time[idx_partner[inter]]-1))))
      }
    }
    behav_freq_table[b,'duration.partner'] = sum(all_logs_partner$duration.s[idx_partner])
    behav_freq_table[b,'proportion.partner'] = behav_freq_table[b,'duration.partner']/session_length[s]
    behav_freq_table[b,'num.events.partner'] = length(idx_partner)
    
    idx_subject = which(all_logs_subject$Behavior == behavs[b] & all_logs_subject$session_num==s)
    sec_subject = vector()
    if(length(idx_subject)>0){
      for (inter in 1:length(idx_subject)){
        
        sec_subject = c(sec_subject, round(c(all_logs_subject$start.time[idx_subject[inter]]:
                                               (all_logs_subject$end.time[idx_subject[inter]]-1))))
      }
    }
    behav_freq_table[b,'duration.subject'] = sum(all_logs_subject$duration.s[idx_subject])
    behav_freq_table[b,'proportion.subject'] = behav_freq_table[b,'duration.subject']/session_length[s]
    behav_freq_table[b,'num.events.subject'] = length(idx_subject)
    
    behav_freq_table[b,'duration.overlap'] = length(intersect(sec_partner, sec_subject))
  }
  
  behav_freq_table_all =rbind(behav_freq_table_all, behav_freq_table)
}

data = behav_freq_table_all[which(behav_freq_table_all$Behavior!= "Proximity"),]

data2 = data.frame()
data2[1,"social_categ"] = "social"
data2[1,"duration"] = sum(data$duration.subject[data$social_categ=="social"]) 
data2[2,"social_categ"] = "asocial"
data2[2,"duration"] = sum(data$duration.subject[data$social_categ=="asocial"]) 
slices <- data2$duration
lbls <- data2$social_categ
pie(slices, labels = lbls, main="Behaviors observed (duration)", cex=2, col=cm.colors(2))
data2$duration[1]/sum(data2$duration)

data3 = data.frame()
data3[1,"valence"] = "affiliative"
data3[1,"duration"] = sum(data$duration.subject[data$valence_categ=="affiliative"]) 
data3[2,"valence"] = "agonistic"
data3[2,"duration"] = sum(data$duration.subject[data$valence_categ=="agonistic"]) 
slices <- data3$duration
lbls <- data3$valence
pie(slices, labels = lbls, main="Behaviors observed (duration)", cex=2, col=heat.colors(2))
data3$duration[1]/sum(data3$duration)

data4 = data.frame()
data4[1,"recip"] = "reciprocal"
data4[1,"duration"] = sum(data$duration.subject[data$recip_categ=="reciprocal"]) 
data4[2,"recip"] = "Not"
data4[2,"duration"] = sum(data$duration.subject[data$recip_categ=="Not"]) 
slices <- data4$duration
lbls <- data4$recip
pie(slices, labels = lbls, main="Behaviors observed (duration)", cex=2, col=terrain.colors(2))
data4$duration[1]/sum(data4$duration)

data5 = data.frame()
data5[1,"overlap"] = "Subject-partner overlap"
data5[1,"duration"] = sum(behav_freq_table_all$duration.overlap)
data5[2,"overlap"] = "Different"
data5[2,"duration"] = sum(session_length) - sum(behav_freq_table_all$duration.overlap)
slices <- data5$duration
lbls <- data5$overlap
pie(slices, labels = lbls, main="Behaviors observed (duration)", cex=2, col=topo.colors(2))
data5$duration[1]/sum(data5$duration)
