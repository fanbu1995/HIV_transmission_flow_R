# 05/23/2021
# make plots for the HIV tranmission flow project (simulation experiments)

library(tidyverse)
setwd("~/Documents/Research_and_References/HIV_transmission_flow/")

dat_weight = read_csv('weight_means.csv')
dat_prop = read_csv('prop_means.csv')

dat_weight$type = 'male source'
dat_prop$type = 'transmission direction'

dat_all = rbind(dat_weight, dat_prop)

#cols = c("#F8766D", "#00BFC4", "#F8766D", "#00BFC4")

## get scenario 1
dat1 = dat_all %>% filter(scenario == 1)

ggplot(data=dat1) +
  geom_violin(aes(y=means, x = type, fill=weight), 
              draw_quantiles = 0.5, show.legend = FALSE) +
  labs(x=NULL, y='post. means', fill='proportions\n of...') +
  #scale_fill_manual(values = cols) + 
  theme_bw(base_size = 14)

## scenario 2
dat2 = dat_all %>% filter(scenario == 2)

ggplot(data=dat2) +
  geom_violin(aes(y=means, x = type, fill=weight), 
              draw_quantiles = 0.5) +
  labs(x=NULL, y='post. means', fill='proportions\n of...') +
  #scale_fill_manual(values = cols) + 
  theme_bw(base_size = 14)


## weights from younger vs older men

### filter out extreme weights
dat_weight = dat_weight %>% 
  filter(means < 0.95 & means > 0.05)

ggplot(data=dat_weight) +
  geom_hline(yintercept = c(0.3,0.6), col = 'brown', 
             size=0.8, linetype = 2) + 
  geom_violin(aes(y=means, x = as.factor(scenario), fill=weight),
              draw_quantiles = 0.5, alpha = 0.7) +
  #geom_boxplot(aes(y=means, x = as.factor(scenario), fill=weight))+
  scale_y_continuous(limits = c(0.1, 0.9)) +
  scale_x_discrete(labels = c('Scenario 1', 'Scenario 2')) + 
  labs(x=NULL, y = 'post. means', fill='male\nsource') +
  theme_bw(base_size = 14)

## props of FM and MF trans
ggplot(data=dat_prop) +
  geom_hline(yintercept = c(0.4, 0.5, 0.6), col = 'brown', 
             size=0.8, linetype = 2) + 
  geom_violin(aes(y=means, x = as.factor(scenario), fill=weight),
              draw_quantiles = 0.5, alpha = 0.7) +
  #geom_boxplot(aes(y=means, x = as.factor(scenario), fill=weight))+
  scale_y_continuous(limits = c(0.3, 0.7)) +
  scale_x_discrete(labels = c('Scenario 1', 'Scenario 2')) + 
  labs(x=NULL, y = 'post. means', fill='transmission\ndirection') +
  theme_bw(base_size = 14)

