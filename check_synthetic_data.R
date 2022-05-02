# 04/24/2020
# check out Olli's synthetic dataset

library(tidyverse)

setwd("~/Downloads/")
dat = read.csv("200423_synthetic data_1.csv",as.is = T)

str(dat)

total_score = dat$PHYLOSCANNER_SCORE_LINKED + 
  dat$PHYLOSCANNER_SCORE_MF + dat$PHYLOSCANNER_SCORE_FM
summary(total_score) # sum: from 1.6xx to 2

summary(dat$PHYLOSCANNER_SCORE_LINKED) # lower-bounded by 0.6
summary(dat$PHYLOSCANNER_SCORE_MF) # from 0 to 1
summary(dat$PHYLOSCANNER_SCORE_FM) # from 0 to 1

hist(dat$PHYLOSCANNER_SCORE_LINKED, breaks = 20)
# The LINKED scores: a uniform(0.6,1) + point mass at 1?

hist(dat$PHYLOSCANNER_SCORE_FM, breaks = 20) 
hist(dat$PHYLOSCANNER_SCORE_MF, breaks = 20)
# The MF, FM scores: a bi-modal distribution??

# 05/01/2020: FM + MF == 1?
summary(dat$PHYLOSCANNER_SCORE_FM + dat$PHYLOSCANNER_SCORE_MF) # all 1!

summary(dat$PHYLOSCANNER_SCORE_MF/dat$PHYLOSCANNER_SCORE_LINKED)
# from 0 to 1.4??? how can it exceed 1??

summary(dat$PHYLOSCANNER_SCORE_FM/dat$PHYLOSCANNER_SCORE_LINKED)
# from 0 to 1.5???

# Question 1: how were the scores generated? (some kind of distribution?)
# Question 1.5: how can MF/FM score exceed LINKED score?
#               how would those be calculated from real phyloscanner data?
#               (cf: lancet-supplement, but vague)