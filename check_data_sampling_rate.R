# 08/26/2021
# take look at real data and sampling rates from Xiaoyue

# 08/31/2021
# try combining estimated surfaces with sampling rates
# (using inland population rates avg. between source and recipient for now...)

# 09/15/2021
# check out unprocessed real data again...
# to see dist. of L and D scores without 0/1 extreme values

library(tidyverse)
library(stringr)

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

# real data (shared by Oliver in June)
setwd('~/Documents/Research_and_References/HIV_transmission_flow/')

dat = read_csv('Rakai_data.csv')
dim(dat)
glimpse(dat)

# check linked score
range(dat$POSTERIOR_SCORE_LINKED)
hist(dat$POSTERIOR_SCORE_LINKED)

# check direction score
range(dat$NETWORK_SCORE_MF)
hist(dat$NETWORK_SCORE_MF)


# filter the data to get "unextreme" L scores and D scores
dat.proc = dat %>% filter(POSTERIOR_SCORE_LINKED > 0.2, 
                          POSTERIOR_SCORE_LINKED < 1,
                          POSTERIOR_SCORE_MF > 0, 
                          POSTERIOR_SCORE_MF < 1)
dim(dat.proc) # 406 rows with scores strictly between 0 and 1
range(dat.proc$POSTERIOR_SCORE_LINKED) #0.2082297 0.9964466
range(dat.proc$POSTERIOR_SCORE_MF) #0.02977233 0.97500000

#get the logit transformed scores
L = logit(dat.proc$POSTERIOR_SCORE_LINKED)
D = logit(dat.proc$POSTERIOR_SCORE_MF)

# try fit GMM on the Ls and Ds
library(mixtools)
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
## GMM for L
GMM.L = normalmixEM2comp(L, lambda = c(0.3, 0.7), mu = c(0, 1),
                         sigsqrd = c(1.5, 1.5))
data.frame(x = GMM.L$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 0.2, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(GMM.L$mu[1], GMM.L$sigma[1], lam = GMM.L$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(GMM.L$mu[2], GMM.L$sigma[2], lam = GMM.L$lambda[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density")

## GMM for D
GMM.D = normalmixEM(D, k = 3)
data.frame(x = GMM.D$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 0.2, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(GMM.D$mu[1], GMM.D$sigma[1], lam = GMM.D$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(GMM.D$mu[2], GMM.D$sigma[2], lam = GMM.D$lambda[2]),
                colour = "blue", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(GMM.D$mu[3], GMM.D$sigma[3], lam = GMM.D$lambda[3]),
                colour = "green", lwd = 1.5) +
  ylab("Density")

# look at the sampling rates shared by Xiaoyue
samp_rates = readRDS('sampling_rate_samples.rds')
glimpse(samp_rates)

unique(samp_rates$SAMPLING_CATEGORY)

# 08/31/2021
# take posterior average of sampling rates within each 
# "SAMPLING_CATEGORY"
avg_samp_rates = samp_rates %>% group_by(SAMPLING_CATEGORY) %>%
  summarise(avg_P = mean(P), avg_LP = mean(LP)) %>%
  select(age_group = SAMPLING_CATEGORY, avg_P, avg_LP)

# get location (inland v.s. fishing), gender, 
#  and age group (take the lower bound only)
avg_samp_rates = avg_samp_rates %>% 
  mutate(location = if_else(str_starts(age_group, 'i:'), 'inland', 'fishing'),
         gender = if_else(str_detect(age_group, ':F:'), 'F', 'M'),
         age = sapply(age_group, 
                      function(s) as.numeric(unlist(str_split(s, '-'))[2])-1)) %>%
  select(location, gender, age, avg_P, avg_LP)

# save it for later use
# write_csv(avg_samp_rates, 'avg_sampling_rates.csv')


# get the cross-product of sampling rates (for each age-age tile)
# use the inland rates for now
F_in = avg_samp_rates %>% filter(gender == 'F', location == 'inland') %>% 
  select(avg_P) %>% pull()
M_in = avg_samp_rates %>% filter(gender == 'M', location == 'inland') %>% 
  select(avg_P) %>% pull()

# outer product to get a sampling rate matrix for each 1-year age tile
cross_rates = outer(M_in, F_in, "*")

# get it to the "wide format" (as a matrix)
cross_rates_dat = as.data.frame(cross_rates)
names(cross_rates_dat) = as.character(c(15:49))
cross_rates_dat$MAge = c(15:49)

# transform it to the "long format"
cross_rates_dat_long = gather(cross_rates_dat, key = 'FAge', value = 'P', -MAge)



###
# a function to invert to the "population" transmission pattern
# and plot it on a surface
inv_trans_surf <- function(dens_file, plot='gg', age_range = c(15:49)+0.5){
  # dens_file: path to the density (eval at midpoints) matrix from Python
  # plot: "gg" for ggplot tile, "base" for base R plot (using image),
  #       FALSE for no plot
  Z =  read_table(dens_file, col_names = FALSE) %>% 
    as.matrix() %>% t() # seems to need transpose??
  Z_adj = Z/cross_rates
  
  # base R plot if...
  if(plot == 'base'){
    
    image(age_range, age_range, Z_adj, xlab='Male age', ylab='Female age')
  }
  
  # transform it to long format
  Z_dat = as.data.frame(Z_adj)
  names(Z_dat) = as.character(age_range)
  Z_dat$MAge = age_range
  Z_dat_long = gather(Z_dat, key = 'FAge', value = 'density', -MAge)
  Z_dat_long$FAge = as.numeric(Z_dat_long$FAge)
  
  # ggplot2 if... 
  if(plot=='gg'){
    print(
    ggplot(Z_dat_long, aes(x=MAge, y=FAge)) + 
      geom_tile(aes(fill = density)) +
      scale_fill_distiller(palette = "Blues", direction = 1) +
      scale_y_continuous(breaks = c(20,30,40,50),
                         expand=c(0,0)) +
      scale_x_continuous(expand=c(0,0)) +
      labs(x = "Male age", y = 'Female age') +
      theme_bw(base_size = 14) +
      theme(legend.position = 'none')
    )
  }
  
  # finally return long formart
  Z_dat_long
}


# do it
# 09/20/2021: do another one
pdf('adj_surfaces_fixAlloc2.pdf', width = 6, height = 6)
## MAP estimates
Z_MF_adj = inv_trans_surf('MF_surface_midpoints_MAP_fixAlloc2.txt')
Z_FM_adj = inv_trans_surf('FM_surface_midpoints_MAP_fixAlloc2.txt')

## Mean surface
Z_MF_mean_adj = inv_trans_surf('MF_surface_midpoints_mean_fixAlloc2.txt')
Z_FM_mean_adj = inv_trans_surf('FM_surface_midpoints_mean_fixAlloc2.txt')

dev.off()

## 09/16/2021
## 09/20/2021: another one
pdf('adj_surfaces_specTreat2.pdf', width = 6, height = 6)
# do this for a version of semi-fixed results (spec. Treat of 0/1 scores)
Z_MF_adj = inv_trans_surf('MF_surface_midpoints_MAP_specTreat2.txt')
Z_FM_adj = inv_trans_surf('FM_surface_midpoints_MAP_specTreat2.txt')

## Mean surface
Z_MF_mean_adj = inv_trans_surf('MF_surface_midpoints_mean_specTreat2.txt')
Z_FM_mean_adj = inv_trans_surf('FM_surface_midpoints_mean_specTreat2.txt')

dev.off()

### plot relative transmission rates for each age group 
#   (for the recipient gender)
#   as lines over the source age

## a helper function to get age bin for each specific age
get_age_bin <- function(x, breaks = seq(15,50,by=1)){
  num_bin = length(breaks) - 1
  res = character(length(x))
  for(i in 1:num_bin){
    lb = breaks[i]
    ub = breaks[i+1]
    res[x>lb & x<ub] = paste0(as.character(lb),'-',as.character(ub))
  }
  res
}

# a function to plot the lines
# 09/20/2021 update: fix the y-axis limits for better comparison
plot_source_freq <- function(Z, recipient = 'F',
                             age_breaks = seq(15,50,by=1), 
                             normalize = TRUE,
                             ylimits = c(0,0.1)){
  # group age first
  Z = Z %>% mutate(MAge_group = get_age_bin(MAge, breaks = age_breaks),
                   FAge_group = get_age_bin(FAge, breaks = age_breaks))
  
  if(recipient=='F'){
    # take the mean density within each recipient age group
    Z = Z %>% group_by(FAge_group, MAge) %>% 
      summarise(mean_dens = mean(density)) %>%
      select(rec_group = FAge_group, source_age = MAge, 
             density=mean_dens)
    
    xlab = 'Male source age'
    colorlab = 'Female Age'
  }else{
    # take the mean density within each recipient age group
    Z = Z %>% group_by(MAge_group, FAge) %>% 
      summarise(mean_dens = mean(density)) %>%
      select(rec_group = MAge_group, source_age = FAge, 
             density=mean_dens)
    
    xlab = 'Female source age'
    colorlab = 'Male Age'
  }
  
  if(normalize){
    # renormalize among sources within each recipient age group
    Z = Z %>% group_by(rec_group) %>%
      mutate(norm_density = density/sum(density)) %>% 
      ungroup() %>% 
      select(rec_group, source_age, density=norm_density)
  }
    
  print(
    ggplot(Z, aes(x=source_age, y=density)) +
      geom_line(aes(color = rec_group)) +
      labs(y='Transmission frequency', x=xlab, color=colorlab) +
      scale_y_continuous(limits = ylimits) +
      theme_bw(base_size = 14) +
      theme(legend.position = 'bottom', 
            legend.title = element_text(size=8),
            legend.text = element_text(size=6),
            legend.key.width = unit(0.3, 'cm'),
            legend.key.height = unit(0.1,'cm')) +
      guides(color = guide_legend(nrow=(length(age_breaks)-1) %/% 12 + 1,
                                  byrow=TRUE))
  )
  
  # return the age grouped data
  #Z
}


# save it to PDF
# order: MAP MF, MAP FM, mean MF, mean FM
# 09/20/2021: do another one

pdf('recipient_source_freq_fixedAlloc2.pdf',
    height=5, width=7)

ylims = c(0,0.11)

plot_source_freq(Z_MF_adj, recipient = 'F', ylimits = ylims)
plot_source_freq(Z_FM_adj, recipient = 'M', ylimits = ylims)
plot_source_freq(Z_MF_mean_adj, recipient = 'F', ylimits = ylims)
plot_source_freq(Z_FM_mean_adj, recipient = 'M', ylimits = ylims)

dev.off()

## 09/16/2021
# do this for the semi-fixed with spec. Treat
## 09/20/2021: another one
pdf('recipient_source_freq_specTreat2.pdf',
    height=5, width=7)

ylims = c(0,0.11)

plot_source_freq(Z_MF_adj, recipient = 'F', ylimits = ylims)
plot_source_freq(Z_FM_adj, recipient = 'M', ylimits = ylims)
plot_source_freq(Z_MF_mean_adj, recipient = 'F', ylimits = ylims)
plot_source_freq(Z_FM_mean_adj, recipient = 'M', ylimits = ylims)

dev.off()

## try 5-year age band
pdf('recipient_source_freq_fixedAlloc_5year.pdf',
    height=5, width=7)

plot_source_freq(Z_MF_adj, recipient = 'F', age_breaks = seq(15,50,by=5))
plot_source_freq(Z_FM_adj, recipient = 'M',age_breaks = seq(15,50,by=5))
plot_source_freq(Z_MF_mean_adj, recipient = 'F',age_breaks = seq(15,50,by=5))
plot_source_freq(Z_FM_mean_adj, recipient = 'M',age_breaks = seq(15,50,by=5))

dev.off()
