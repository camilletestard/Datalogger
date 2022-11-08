library(gganimate)

setwd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/Amos_2021-07-29'); 

umap_data = read.csv2('umap_output.csv',sep=',',header=FALSE)
colnames(umap_data)=c("dim1","dim2","label","window")

# query.points = data.frame(
# 	i = 'stage',
# 	x = x,
# 	y = y
# )
# 
# 
# p = ggplot() + 
# 	geom_point(data=ref.points,aes(x,y)) +
# 	geom_point(data=query.points,aes(x,y)) +


p = plot.umap(umaps.all,color=color.by,file=NULL,legend=FALSE) +
  ggtitle('{closest_state}') +
  transition_states(iter, transition_length = 4, state_length = 1) +
  ease_aes('cubic-in-out')
push.status('plot.umap')

a = animate(
  plot = p,
  #	renderer = gifski_renderer(loop = FALSE),
  nframes = 200, width = 800, height = 640, res = 300
)
push.status('animate')

anim_save(
  file.path('animations',paste0('umap-harmony-',prefix,'-',gsub('_','',color.by),'-all_nolegend.gif')),
  a
)