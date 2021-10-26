library(tidyverse)
library(tidyquant)
library(reshape2)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)


data <- read_csv('../data/zach_results.csv')
data <- data %>% mutate(flag = ifelse(p_l == 0 | p_l == 1, 'one', 'two'))
dev.off()
p1 = ggplot(data = data,
       aes(x = t, y = p_ab,
           color = as.factor(p_l)))+
  geom_line(alpha = 1)+
  theme(aspect.ratio = 3/4,
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            legend.position = 'none',
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            #axis.text = element_blank(),
            axis.ticks = element_blank())+
  scale_y_log10(limits = c(1e-39, 1))+
  scale_color_brewer(palette = 'Set1', direction = -1)+
  labs(y = "Abundance",
       x = "Time")
print(p1)

p2 = ggplot(data = data,
            aes(x = t, y = p_ab_per))+
  geom_line(aes(linetype = flag,
                color = as.factor(p_l)))+
  theme(aspect.ratio = 3/4,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = 'none',
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_linetype_manual(values = c('dashed', 'solid'))+
  scale_y_log10(limits = c(1e-39, 1))+
  scale_color_brewer(palette = 'Set1', direction = -1)+
  labs(y = "",
       x = "Time")
print(p2)

#Load matrices
load('../data/example_perturbation_parameters.RData')

A_pre = my_list$A_pre
A_pre_melt = melt(A_pre)

A_pre_plot =  ggplot(data = A_pre_melt, 
                     aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile()+
  theme(panel.background = element_rect(fill = NA),
        plot.background = element_blank(),
        # legend.title = element_text(size = 15, hjust = 0.2,
        #                             margin = margin(b = -10)),
        legend.title = element_blank(),
        legend.position = c(1.16, 0.55),
        legend.background = element_rect(fill = NA),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 1)+
  labs(x = '', y = '',
       fill = expression(A[ij]))+
  scale_fill_distiller(palette = 'RdYlBu',
                       guide = guide_colorbar(barwidth = 0.6,
                                              barheight = 4.5),
                       breaks = c(0.5, 0.9, 1.3))
print(A_pre_plot)

A_post = my_list$A_post
A_post_melt = melt(A_post)

A_post_plot =  ggplot(data = A_post_melt, 
                     aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile()+
  theme(panel.background = element_rect(fill = NA),
        plot.background = element_blank(),
        legend.position = c(1.16, 0.55),
        legend.background = element_rect(fill = NA),
        # legend.title = element_text(size = 15, hjust = 0.2,
        #                             margin = margin(b = -10)),
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 1)+
  labs(x = '', y = '',
       fill = expression(A[ij]))+
  scale_fill_distiller(palette = 'RdYlBu',
                       guide = guide_colorbar(barwidth = 0.6,
                                              barheight = 4.5),
                       breaks = c(0.5, 0.9, 1.3))
  

print(A_post_plot)


data = data.frame(x = as.vector(A_pre), y = as.vector(A_post))

r_2 = paste("R^2 == ", round(summary(lm(data$y~data$x))$r.squared, digits = 3))

correlation_plot = ggplot(data = data,
                          aes(x = x, y = y))+
  geom_point()+
  geom_abline(slope = 1)+
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        legend.position = 'none',
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text = element_text(size = 10))+
  scale_x_continuous(breaks = c(0.6, 0.9, 1.2, 1.5),
                     limits = c(0.4, 1.5))+
  scale_y_continuous(breaks = c(0.6, 0.9, 1.2, 1.5),
                     limits = c(0.4, 1.5))+
  annotate('text', x = 0.8, y = 1.4, label = r_2, parse = T)+
  labs(x = 'A', y = expression(A[p]))
  
print(correlation_plot)

yheight = 0.63
xheight = -0.30
height = 0.34
ptot1 = ggdraw() +
  draw_plot(p1) +
  draw_plot(A_pre_plot, x = xheight, y = yheight, height = height)

ptot2 = ggdraw() +
  draw_plot(p2) +
  draw_plot(A_post_plot, x = xheight, y = yheight, height = height)

lay <- rbind(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3),
             c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3),
             c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3))
pdf('../sandbox/dynamics_collapse.pdf', height = 3.5, width  = 12)
grid.arrange(ptot1, ptot2, correlation_plot,
             layout_matrix = lay)
dev.off()

plot1_bis = 
ggdraw()+
  draw_plot(p1)+ 
  draw_plot(correlation_plot, x = 0.2, y = 0.32, height = 0.3
            )

lay <- rbind(c(1, 1, 1, 1, 2, 2, 2, 2),
             c(1, 1, 1, 1, 2, 2, 2, 2),
             c(1, 1, 1, 1, 2, 2, 2, 2))
pdf('../sandbox/dynamics_collapse2.pdf', height = 5.53, width  = 8.45)
grid.arrange(plot1_bis, p2,
             layout_matrix = lay)
dev.off()



