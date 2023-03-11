

library(ggplot2)
library(lubridate)

setwd('~/Dropbox (Penn)/Datalogger/Results/All_sessions/Number of units/')
data=read.csv('Session_log_num_units_final.csv')
data$monkey=as.factor(data$monkey)
#data$date=ymd(data$date)

setwd(paste('~/Dropbox (Penn)/Datalogger/Results/'))
ggplot(data = data, aes(x=date, y=num.units, color = monkey))+
  geom_point(size = 3)+ theme_classic(base_size=16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylim(200, 350)+ylab('Number of units')
ggsave('Number of units for all sessions.pdf')

setwd(paste('~/Dropbox (Penn)/Datalogger/Results/'))
ggplot(data = data, aes(x=date, y=num.units.TEO))+
  geom_point(size = 3, color='red')+ 
  geom_point(aes(x=date, y=num.units.vlPFC), size = 3, color='blue')+
  theme_classic(base_size=16)+
  facet_wrap(~monkey)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylim(100, 200)+ylab('Number of units')
ggsave('Number of units split by area.pdf')
