# 08/31/2021
# HIV transmission utils
# some helper functions 
# AND real data plotting
# define logit and expit functions

library(tidyverse)
library(ggplot2)

setwd('~/Documents/Research_and_References/HIV_transmission_flow/')

# logit
logit <- function(p){
  # p should be between 0 and 1
  log(p/(1-p))
}

# expit
expit <- function(x){
  # x is a real number
  exp(x)/(1+exp(x))
}

# real data peek
dat = read_csv('Rakai_data_processed.csv')

# create pdf document for saving these plots
pdf('real_data_plots.pdf', width = 8, height = 4)

# look at linkage score distribution (truncated at L > 0.2)
ggplot(dat) + 
  geom_histogram(aes(x=POSTERIOR_SCORE_LINKED)) +
  labs(x='Linkage score (L)') +
  theme_bw(base_size = 14)

# look at MF direction score
ggplot(dat) +
  geom_histogram(aes(x=POSTERIOR_SCORE_MF), fill='skyblue') +
  labs(x='M->F direction score (D)') +
  theme_bw(base_size = 14)

# MF direction score and linkage scores 
# on the logit transformed scale
pdf('real_data_plots_transform.pdf', 
    width = 8, height = 4)
dat = dat %>% mutate(transD = logit(POSTERIOR_SCORE_MF),
                     transL = logit(POSTERIOR_SCORE_LINKED)) 
ggplot(dat) +
  geom_histogram(aes(x=transD), fill='skyblue') +
  labs(x='M->F direction score (D), transformed') +
  theme_bw(base_size = 14)

ggplot(dat) +
  geom_histogram(aes(x=transL), fill='skyblue') +
  labs(x='Linkage score (L), transformed') +
  theme_bw(base_size = 14)

## (if filter at higher thres. for L??)
# dat1 = dat %>% filter(POSTERIOR_SCORE_LINKED > 0.4)
# ggplot(dat1) +
#   geom_histogram(aes(x=transD), fill='skyblue') +
#   labs(x='M->F direction score (D), transformed') +
#   theme_bw(base_size = 14)

# check for D scores and L scores at the boundaries (0, 1)
dat = dat %>%
  mutate(D_extreme = 
           if_else(POSTERIOR_SCORE_MF == 0.98 | POSTERIOR_SCORE_MF == 0.02,
                   1, 0),
         L_extreme = 
           if_else(POSTERIOR_SCORE_LINKED == 0.98, 1, 0))

# histogram of L scores with extreme D scores
ggplot(dat %>% filter(D_extreme==1)) +
  geom_histogram(aes(x=POSTERIOR_SCORE_LINKED), fill='skyblue', bins=15) +
  labs(x='Linkage score (L), with 0/1 D scores') +
  theme_bw(base_size = 14)

# histogram of D scores with L=1
ggplot(dat %>% filter(L_extreme==1)) +
  geom_histogram(aes(x=POSTERIOR_SCORE_MF), fill='skyblue', bins=15) +
  labs(x='MF direction score (D), with L = 1') +
  theme_bw(base_size = 14)

# Age pair distribution
## get types by fixed thresholds...
## L < 0.6: 0 surface
## L >= 0.6 & D > 0.5: MF
## L >= 0.6 & D <=0.5: FM

dat = dat %>% 
  mutate(type = case_when(POSTERIOR_SCORE_LINKED < 0.6 ~ '0',
                          POSTERIOR_SCORE_LINKED >= 0.6 & 
                            POSTERIOR_SCORE_MF > 0.5 ~ "MF",
                          TRUE ~ 'FM'))

# hist of direction score colored by types
ggplot(dat) +
  geom_histogram(aes(x=dat$POSTERIOR_SCORE_MF, fill = type)) +
  labs(x='M->F direction score (D)') +
  theme_bw(base_size = 14)

# 09/25/2021
# hist of direction scores faceted by types
ggplot(dat) +
  geom_histogram(aes(x=dat$POSTERIOR_SCORE_MF)) +
  labs(x='M->F direction score (D)') +
  facet_grid(rows = vars(type)) + 
  theme_bw(base_size = 14)

# distribution of age pairs colored by types
ggplot(dat) +
  geom_point(aes(x=MALE_AGE_AT_MID, y=FEMALE_AGE_AT_MID, color=type)) +
  labs(x='Male age', y='Female age') +
  theme_bw(base_size = 14)

# distribution of age pairs faceted by types       
ggplot(dat) +
  geom_point(aes(x=MALE_AGE_AT_MID, y=FEMALE_AGE_AT_MID)) +
  labs(x='Male age', y='Female age') +
  facet_grid(cols = vars(type)) +
  theme_bw(base_size = 14)

# slipt pairs by if they have extreme D or L scores
ggplot(dat %>% filter(D_extreme == 1 | L_extreme == 1)) +
  geom_point(aes(x=MALE_AGE_AT_MID, y=FEMALE_AGE_AT_MID)) +
  labs(x='Male age', y='Female age', title = 'with 0/1 scores') +
  facet_grid(cols = vars(type)) +
  theme_bw(base_size = 14)

ggplot(dat %>% filter(D_extreme == 0, L_extreme == 0)) +
  geom_point(aes(x=MALE_AGE_AT_MID, y=FEMALE_AGE_AT_MID)) +
  labs(x='Male age', y='Female age', title = 'with non 0/1 scores') +
  facet_grid(cols = vars(type)) +
  theme_bw(base_size = 14)

dev.off()
