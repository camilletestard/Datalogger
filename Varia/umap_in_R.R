library(tidyverse)
library(ggplot2)
library(umap)
library(RColorBrewer)
library(R.matlab)

#Load data:
file = file.choose()  
data <- readMat(file)
neural.data<- as.data.frame(data$Spike.count.raster)
labels<- as.factor(data$behavior.labels)
hist(table(labels), freq=TRUE, xlab = levels(labels), ylab = "Frequencies")

#Select behaviors of interest
behav = c(4:6,17)
idx = which(labels %in% behav)
neural.data.final = neural.data[idx,]
labels.final = labels[idx]

set.seed(3)

a = umap(neural.data.final)#, n_neighbors = 5, min_dist = 0.5)
a.umap = data.frame(as.data.frame(a$layout),labels.final)

ggplot(a.umap,aes_string('V1','V2',color='labels.final')) +
  geom_point(size=1.25) +
  #scale_color_manual(values=region.colors) +
  theme_classic(base_size=18) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  coord_fixed() +
  theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap_batched.pdf',useDingbats=FALSE)