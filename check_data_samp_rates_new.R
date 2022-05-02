# 1/11/2022
# look at real data again
# (data shared in Nov 2021)

# 1/12/2022
# look at de-bugged data
# (updated in Jan 2022)

library(tidyverse)
library(stringr)
library(zoo)
#require(Hmisc)

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

# real data (shared by Oliver in Nov 2021)
setwd('~/Documents/Research_and_References/HIV_transmission_flow/')

#dat = read_csv('211126_data_not_unlike_real_data.csv')
dat = read_csv('Rakai_data_Jan2022.csv')
dim(dat)
glimpse(dat)

# the sampling rates (by individuals)
samp_rates = read_csv('211126_sampling_data.csv')
glimpse(samp_rates)
### for each person and each category (recipient/source)
### there are 500 monte carlo samples

samp_rates = samp_rates %>% 
  mutate(CATE = if_else(str_starts(TYPE, 'REC'), 'REC', 'SC')) %>%
  select(AID = AID2, CATE, MONTE_CARLO_ID, S_PROB)
# CATE: REC for recipient, SC for source


# cross ref. with the data to get age for each person
## get a data frame with AID and age only
id_age = data.frame(AID = c(dat$FEMALE_AID, dat$MALE_AID),
                    Age = c(dat$FEMALE_AGE_AT_MID, dat$MALE_AGE_AT_MID)) %>%
  distinct()

samp_rates_age = left_join(samp_rates, id_age, by = "AID")

## sanity check to see if same num of rows are preserved
nrow(samp_rates)
nrow(samp_rates_age)

# Jan 12 2022: everything fixed now!!

# ## check what's wrong...
# for(id in unique(samp_rates$AID)){
#   rows_id = samp_rates_age %>% filter(AID == id) %>% nrow()
#   
#   if(rows_id != 1000){
#     cat('something wrong with id', id, '!\n')
#   }
# }
# 
# ## turns out some people have more than one different ages!
# ## in fact, some have two or three unique age values
# ## but those ages are close...
# for(id in unique(id_age$AID)){
#   ages = id_age %>% filter(AID == id) %>% select(Age) %>% pull()
#   if(length(ages) > 1){
#     cat('Multiple ages for id', id, ":", paste(ages, collapse = ','),'.\n')
#   }
# }
# 
# ## what I could do:
# ## randomly assign one unique age value to each person
# ## otherwise the whole thing breaks down
# 
# ## BUT to keep it simple (since I'm mapping it back to points on age space anyway)
# ## just work with the multiple ages for some people

## save it for now
saveRDS(samp_rates_age, './samp_rates_age.rds')


# 01/12/2022:
# get the monte carlo samples of sampling rates within each 1-yr age band
samp_rates_age = readRDS('./samp_rates_age.rds') %>%
  mutate(gender = if_else(str_starts(AID, 'F'), 'F', 'M'),
         age_group = floor(Age)) # age_group: lower bound of 1-yr age band

## save this one instead
saveRDS(samp_rates_age, './samp_rates_age_gender.rds')

## a function to randomly select a certain number of samples for
## age bands of a certain gender for a category
get_samps_sel <- function(rates, num_samp = 200, 
                          age_range = c(15, 49),
                          gen = 'M', cate = 'SC',
                          interpolate = TRUE){
  rates_sel = rates %>% 
    filter(gender == gen, CATE == cate, 
           (age_group >= age_range[1]) & (age_group <= age_range[2])) %>%
    select(gender, CATE, prob = S_PROB, age_group) %>%
    group_by(age_group) %>%
    sample_n(num_samp) %>%
    ungroup() %>%
    arrange(age_group)
  
  n_age = length(unique(rates_sel$age_group))
  rates_sel$monte_carlo_id = rep(1:num_samp, n_age)
  
  # use interpolation to fill in the gaps
  if(interpolate){
    #N = (age_range[2] - age_range[1] + 1) * num_samp
    template = data.frame(gender = gen, CATE = cate, 
                          age_group = rep(age_range[1]:age_range[2], num_samp)) %>%
      arrange(age_group)
    n_age = age_range[2] - age_range[1] + 1
    template$monte_carlo_id = rep(1:num_samp, n_age)
    res = left_join(template, rates_sel) %>%
      group_by(monte_carlo_id) %>%
      #mutate(approx_prob = na.approx(prob)) %>% # linear interpolation doesn't work!
      #mutate(approx_prob = na.spline(prob)) %>% # splines can to extrap. but values seem a bit extreme
      mutate(approx_prob = Hmisc::approxExtrap(x = age_group,
                                               y = prob,
                                               xout = age_group)$y) %>% # try this linear extrap. thing
      ungroup() %>%
      select(gender, CATE, age_group, monte_carlo_id, prob = approx_prob)
  }else{
    # n_age = length(unique(rates_sel$age_group))
    # rates_sel$monte_carlo_id = rep(1:num_samp, n_age)
    res = rates_sel %>% select(gender, CATE, age_group, monte_carlo_id, prob)
  }
  
  res
}

## try to generate a male source sampling rates 
M_SC_samps = get_samps_sel(samp_rates_age, num_samp = 200)
F_SC_samps = get_samps_sel(samp_rates_age, num_samp = 200, gen = 'F')

## save those source sampling rates
saveRDS(list(M_SC_samps = M_SC_samps, F_SC_samps = F_SC_samps),
        'source_samp_rates_age.rds')

# 01/11/2022 stuff
# Simple handling of the sampling rates
# Averaging over the monte carlo samples 
# AND create new field for gender
samp_rates_age = samp_rates_age %>% group_by(AID, CATE, Age) %>%
  summarise(Avg_prob = mean(S_PROB)) %>%
  ungroup()

samp_rates_age = samp_rates_age %>%
  mutate(gender = if_else(str_starts(AID, 'F'), 'F', 'M'),
         age_group = floor(Age)) # age_group: lower bound of 1-yr age band
samp_rates_age %>% filter(gender == 'F') %>% select(age_group) %>% 
  pull() %>% unique() %>%sort()
samp_rates_age %>% filter(gender == 'M') %>% select(age_group) %>% 
  pull() %>% unique() %>%sort()
## female age group: 14, then 16-49 (an extra 14 something)
## male age group: 18-51 (extra 51 something, and missing 15, 16 & 48 groups)

## average by gender, category and age groups
samp_rates_group = samp_rates_age %>% group_by(CATE, gender, age_group) %>%
  summarise(mean_prob = mean(Avg_prob)) %>%
  ungroup()

# some interpolation to fill in the age gaps 
# (making sure each age band within the range has a value)
## a function to do that (via splines for now)
## (o.w. can use na.approx for linear interpolation)

nice_interpolate <- function(rates){
  age_range = range(rates$age_group)
  age_bands = c(age_range[1]:age_range[2])
  n = length(age_bands)
  template = data.frame(CATE = rep(c('REC', 'SC'), each = n*2),
                        gender = rep(c('F','M','F','M'), each = n),
                        age_group = rep(age_bands, 4))
  res = left_join(template, rates) %>%
    mutate(approx_prob = na.spline(mean_prob)) %>%
    select(CATE, gender, age_group, prob = approx_prob)
  res
}

samp_rates_smooth = nice_interpolate(samp_rates_group)

# then produce the cross product sampling rates for each age tile
# stratified by M->F and F->M surface
# output: 
get_crosstab <- function(rates){
  rl = list()
  ages = c(range(rates$age_group)[1]:range(rates$age_group)[2])
  
  ## MF surface
  F_rec = rates %>% filter(CATE == 'REC', gender == 'F') %>%
    select(prob) %>% pull()
  M_sc = rates %>% filter(CATE == 'SC', gender == 'M') %>%
    select(prob) %>% pull()
  MF_cross = outer(M_sc, F_rec, "*") # M in the rows, F in the cols
  rownames(MF_cross) = colnames(MF_cross) = as.character(ages)
  rl[['MF']] = MF_cross
  
  ## FM surface
  M_rec = rates %>% filter(CATE == 'REC', gender == 'M') %>%
    select(prob) %>% pull()
  F_sc = rates %>% filter(CATE == 'SC', gender == 'F') %>%
    select(prob) %>% pull()
  FM_cross = outer(F_sc, M_rec, "*") # F in the rows, M in the cols
  rownames(FM_cross) = colnames(FM_cross) = as.character(ages)
  rl[['FM']] = FM_cross
  
  rl
}

cross_rates = get_crosstab(samp_rates_smooth)

## save the cross rates
saveRDS(cross_rates, 'cross_samp_rates_new.rds')
