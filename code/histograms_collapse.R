library(tidyverse)
library(scale)
library(gridExtra)
library(cowplot)
setwd("~/Desktop/phd/plant_soil_feedback/code/")
data = read.csv('../data/results.csv', row.names = 1)
dev.off()
#Remove non-convergent data
data = data %>% filter(n_p != -1)
#Get number of different initial species
n_init_vec = unique(data$n_p)
#Get total number of simulations
n_sim = max(data$n_sim)
#Initialize list of plots
myplots = list()
#Plot histograms for each value of n_init_vec
for (i in seq(length(n_init_vec))){
  #Prepare data
  data_i = data %>% 
    filter(n_p == n_init_vec[i]) %>% 
    gather(key = 'n_sp', value = 'count', -n_sim)
  #Generate plot
  p_i = ggplot(data = data_i, aes(x = count, 
                                  fill = n_sp, 
                                  color = n_sp,
                                  linetype = n_sp))+ 
    geom_histogram(position = 'identity', 
                   alpha = 0.5,
                   binwidth = 1)+
    #coord_cartesian(xlim = c(1.5, 6.5))+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 15),
          legend.position = 'top',
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12, hjust = 0.1),
          legend.title.align = 0.5,
          legend.box.margin = margin(b = -50),
          panel.background = element_blank(), 
          panel.border = element_rect(color = 'black', fill=NA),
          aspect.ratio = 0.5)+
    scale_fill_manual(values = c('red', 'blue'),
                      labels = c('Before', 'After'),
                      guide = guide_legend(direction = "horizontal",
                                           title.position = "top"))+
    scale_color_manual(values = c('red', 'blue'),
                       labels = c('Before', 'After'))+
    scale_linetype_manual(values = c('solid', 'dashed'),
                          labels = c('Before', 'After'))+
    scale_x_continuous(breaks = n_init_vec,
                       limits = c(1.5, 6.5))+
    scale_y_continuous(breaks = seq(0, n_sim + 1, length.out = n_sim/3),
                       limits = c(0, n_sim + 1),
                       expand = expansion(mult = c(0, .1)))+
    labs(fill = 'Number of coexisting species',
         color = 'Number of coexisting species',
         linetype = 'Number of coexisting species')
  legend <- get_legend(p_i)
  #Store
  myplots[[i]] =  p_i + theme(legend.position = "none")
}

grid.arrange(legend, myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], nrow = length(n_init_vec) + 1)
ggsave('../sandbox/histogram_collapse.pdf')

