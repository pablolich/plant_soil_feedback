library(tidyverse)

data = read_csv('../data/results.csv')
n_sim = round(max(data$n_sim), digits = -2) # number of simulations per inital richness

# Remove non-convergent data (need to filter on n_p_f, rather than n_p)
# Hopefully we can deal with these cases in the simulation code and eventually
# remove this line
#data = data %>% filter(n_p_f != -1)

# Reshape data for easier ggplotting
data = data %>% 
  mutate(n_p_s = n_p) %>%
  pivot_longer(cols = c(n_p_s, n_p_f),
               names_to = "type",
               values_to = "num_species")

data$type = factor(data$type, levels = c('n_p_s', 'n_p_f'))

data %>% ggplot(aes(x = num_species,
                    after_stat(density),
                    fill = type,
                    color = type,
                    linetype = type)) + 
  geom_histogram(position = "identity",
                 alpha = 0.5,
                 binwidth = 1,
                 size = 0.3) +
  facet_wrap(.~n_p, nrow = 4) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 15),
         strip.text = element_text(size = 15),
         legend.title = element_blank(),
         legend.position = c(0.78, 0.94),
         legend.text = element_text(size = 12, hjust = 0.1),
         panel.background = element_blank(), 
         panel.border = element_rect(color = 'black', fill = NA), 
         aspect.ratio = 0.5) + 
  scale_fill_manual(values = c('blue', 'red'),
                    labels = c('Initial', 'Final')) + 
  scale_color_manual(values = c('blue', 'red'),
                     labels = c('Initial', 'Final'))+
  scale_linetype_manual(values = c('dashed', 'solid'),
                        labels = c('Initial', 'Final'))+
  scale_x_continuous(breaks = seq(1:max(data$num_species)),
                     limits = c(0.5, max(data$num_species) + 0.5)) +
  scale_y_continuous(breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0, 0.1))) + 
  xlab("Number of species") + 
  ylab("Proportion of communities")

ggsave('../sandbox/histogram_collapse_feas.pdf')
