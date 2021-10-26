library(tidyverse)
library(tidyquant)
library(reshape2)
library(cowplot)
library(gridExtra)


data = read_csv('../data/zach_results.csv')
dev.off()
p1 = ggplot(data = data,
       aes(x = t, y = p_ab,
           color = as.factor(p_l)))+
  geom_line(alpha = 0.2)+
  geom_ma(ma_fun = EMA, n = 150, linetype = 'solid')+
theme(aspect.ratio = 3/4,
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            legend.position = 'none',
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.text = element_blank(),
            axis.ticks = element_blank())+
  labs(y = "Abundance",
       x = "Time")
print(p1)

p2 = ggplot(data = data,
            aes(x = t, y = p_ab_per,
                color = as.factor(p_l)))+
  geom_line(alpha = 0.2)+
  geom_ma(ma_fun = EMA, n = 150, linetype = 'solid')+
  theme(aspect.ratio = 3/4,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = 'none',
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
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
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 1)+
    labs(x = '', y = '')

print(A_pre_plot)

A_post = my_list$A_post
A_post_melt = melt(A_post)

A_post_plot =  ggplot(data = A_post_melt, 
                     aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile()+
  theme(panel.background = element_rect(fill = NA),
        plot.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 1)+
  labs(x = '', y = '')

print(A_post_plot)


data = data.frame(x = as.vector(A_pre), y = as.vector(A_post))

correlation_plot = ggplot(data = data,
                          aes(x = x, y = y))+
  geom_point()+
  geom_abline(slope = 1)+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = 'none',
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text = element_text(size = 15))+
  labs(x = 'A', y = expression(A[p]))
  
print(correlation_plot)

yheight = 0.6
xheight = -0.335
height = 0.15
ptot1 = ggdraw() +
  draw_plot(p1) +
  draw_plot(A_pre_plot, x = xheight, y = yheight, height = height)

ptot2 = ggdraw() +
  draw_plot(p2) +
  draw_plot(A_post_plot, x = xheight, y = yheight, height = height)

lay <- rbind(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3),
             c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3),
             c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3))
grid.arrange(ptot1, ptot2, correlation_plot,
             layout_matrix = lay)


