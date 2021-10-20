
library(tidyverse)
setwd("~/Desktop/phd/plant_soil_feedback/code/")
data = read.csv('../data/results.csv', row.names = 1)
#Plot histogram for 2 species 
#Get data 2 sp
data_2sp = data %>% filter(n_p == 2) %>% gather(key = 'n_sp', value = 'count', -n_sim)
ggplot(data = data_2sp, aes(x = count, fill = n_sp))+ 
  geom_histogram(position = 'dodge')
data_3sp = data %>% filter(n_p == 3) %>% gather(key = 'n_sp', value = 'count', -n_sim)
ggplot(data = data_3sp, aes(x = count, fill = n_sp))+ 
  geom_histogram(position = 'dodge')
data_5sp = data %>% filter(n_p == 5) %>% gather(key = 'n_sp', value = 'count', -n_sim)
ggplot(data = data_5sp, aes(x = count, fill = n_sp))+ 
  geom_histogram(position = 'dodge')
data_6sp = data %>% filter(n_p == 6) %>% gather(key = 'n_sp', value = 'count', -n_sim)
ggplot(data = data_6sp, aes(x = count, fill = n_sp))+ 
  geom_histogram(position = 'dodge')

