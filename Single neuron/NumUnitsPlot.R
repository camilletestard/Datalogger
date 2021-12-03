setwd('~/Dropbox (Penn)/Datalogger/Results/')
data=read.csv('Session_log_num_units.csv')
data$monkey=as.factor(data$monkey)


library(ggplot2)

setwd(paste('~/Dropbox (Penn)/Datalogger/Results/', sessions[s],"/",sep=""))
ggplot(data = data, aes(x=date, y=num.units, color = monkey))+
  geom_point(size = 3)+ theme_classic(base_size=16)+
  theme(axis.text.x = element_text(angle = 90))+
  ylim(200, 350)+ylab('Number of units')
ggsave('Number of units for all sessions.png')