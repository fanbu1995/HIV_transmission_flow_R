# 05/17/2022
# further processing real data 
# trying to resolve "MF and FM scores don't sum to 1" issue

library(tidyverse)

dat = read_csv("Rakai_data_Jan2022.csv")

prob_inds = which(dat$POSTERIOR_SCORE_MF + dat$POSTERIOR_SCORE_FM < 0.99)

dat %>% slice(prob_inds) %>%
  select(POSTERIOR_SCORE_MF, POSTERIOR_SCORE_FM) %>%
  print(n=25)

# 2 things:
# if both MF and FM are 0; throw them away
# if only one is 0, then set the other one to 1

dat = dat %>% 
  filter(POSTERIOR_SCORE_MF + POSTERIOR_SCORE_FM > 0)

boost_MF_inds = which((dat$POSTERIOR_SCORE_MF + dat$POSTERIOR_SCORE_FM < 0.99) &
                        (dat$POSTERIOR_SCORE_FM == 0))

boost_FM_inds = which((dat$POSTERIOR_SCORE_MF + dat$POSTERIOR_SCORE_FM < 0.99) &
                        (dat$POSTERIOR_SCORE_MF == 0))

dat$POSTERIOR_SCORE_MF[boost_MF_inds] = 1
dat$POSTERIOR_SCORE_FM[boost_FM_inds] = 1

# check again
prob_inds = which(dat$POSTERIOR_SCORE_MF + dat$POSTERIOR_SCORE_FM < 0.99)

dat %>% slice(prob_inds) %>%
  select(POSTERIOR_SCORE_MF, POSTERIOR_SCORE_FM)

# -- we are good to go!!

write_csv(dat, file = 'Rakai_data_Jan2022_processed.csv')
