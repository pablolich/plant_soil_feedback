library(tidyverse)
library(tidyquant)
library(reshape2)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)
  
  
  data <- read_csv('../data/zach_results.csv')
  data <- data %>% mutate(flag = ifelse(p_l == 0 | p_l == 1, 'one', 'two'))
  colnames(data) <- c('X1', 'Perturbed', 't', 'p_l', 'Unperturbed', 'flag')
  
  data = data %>% 
    pivot_longer(cols = c(Unperturbed, Perturbed),
                 names_to = "type",
                 values_to = "abundance")
  
  #dev.off()
  p1 = ggplot(data = data,
         aes(x = t, y = abundance,
             color = as.factor(p_l)))+
    geom_line(alpha = 0.8, size = 1)+
    facet_wrap(.~factor(type, levels = rev(levels(factor(type)))))+
    theme_classic() + 
    theme(aspect.ratio = 3/4,
              panel.background = element_blank(),
              panel.border = element_rect(fill = NA),
              legend.position = 'none',
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              strip.text = element_text(size = 12))+
    scale_y_log10(limits = c(10^-6, 1), 
                  breaks = c(10^-1, 10^-3, 10^-5), 
                  labels = c(expression(10^-1), expression(10^-3), expression(10^-5)))+
    scale_color_brewer(palette = 'Set1', direction = -1)+
    labs(y = "Frequency",
         x = "Time")
  print(p1)
  
  #Load matrices
  load('../data/example_perturbation_parameters.RData')
  
  A_pre = my_list$A_pre
  A_post = my_list$A_post
  
  data = data.frame(x = as.vector(A_pre), y = as.vector(A_post))
  
  r_2 = paste("R^2 == ", round(summary(lm(data$y~data$x))$r.squared, digits = 3))
  
  correlation_plot = ggplot(data = data,
                            aes(x = x, y = y))+
    geom_point()+
    geom_abline(slope = 1)+
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = 'transparent',
                                          colour = NA),
          plot.background = element_rect(fill = 'transparent', 
                                         colour = NA),
          panel.border = element_rect(fill = 'transparent'),
          legend.position = 'none',
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 7))+
    scale_x_continuous(breaks = c(0.6, 0.9, 1.2, 1.5),
                       limits = c(0.4, 1.5))+
    scale_y_continuous(breaks = c(0.6, 0.9, 1.2, 1.5),
                       limits = c(0.4, 1.5))+
    annotate('text', x = 0.8, y = 1.4, label = r_2, parse = T, size = 3)+
    labs(x = expression(A[u]), y = expression(A[p]))
  
  plot1_bis = 
  ggdraw()+
    draw_plot(p1)+ 
    draw_plot(correlation_plot, x = -0.07, y = 0.24, height = 0.28)
  
  
  pdf('../sandbox/dynamics_collapse2_test.pdf', height = 5.53, width  = 8.45)
  grid.arrange(plot1_bis)
  dev.off()
  
  
  
