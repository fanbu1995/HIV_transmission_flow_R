# 01/12/2022
# post process real data analysis results

library(tidyverse)

## old Macbook path
setwd('~/Documents/Research_and_References/HIV_transmission_flow/')

## new Macbook path
setwd("/Users/fan/Documents/Research/HIV_transmission_flow")

## read cross rates processed using updated data
cross_rates = readRDS('cross_samp_rates_new.rds')
MF_rates = cross_rates$MF
FM_rates = cross_rates$FM

#####
# copied from previous script
# a function to invert to the "population" transmission pattern
# and plot it on a surface
## 01/15/2022 update: trim the cross rates depending on the age range

inv_trans_surf <- function(dens_file, rates,
                           plot='gg', age_range = c(15:49)+0.5){
  # dens_file: path to the density (eval at midpoints) matrix from Python
  # rates: the cross-tile sampling rates
  # plot: "gg" for ggplot tile, "base" for base R plot (using image),
  #       FALSE for no plot
  Z =  read_table(dens_file, col_names = FALSE) %>% 
    as.matrix() %>% t() # seems to need transpose
  
  # trim the rates
  ages = floor(range(age_range))
  st = which(colnames(rates) == ages[1])
  en = which(colnames(rates) == ages[2])
  rates = rates[st:en,st:en]
  
  # adjust
  Z_adj = Z/rates
  
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
  
  # finally return long format
  Z_dat_long
}

## adjust results
# 1. fixed alloc case
## MAP estimates
Z_MF_adj = inv_trans_surf('MF_surface_midpoints_MAP_fixAlloc_Jan2022.txt', MF_rates)
Z_FM_adj = inv_trans_surf('FM_surface_midpoints_MAP_fixAlloc_Jan2022.txt', FM_rates)

## Mean surface
Z_MF_mean_adj = inv_trans_surf('MF_surface_midpoints_mean_fixAlloc_Jan2022.txt', MF_rates)
Z_FM_mean_adj = inv_trans_surf('FM_surface_midpoints_mean_fixAlloc_Jan2022.txt', FM_rates)

# 2. flexible point labels
Z_MF_adj = inv_trans_surf('MF_surface_midpoints_MAP_specTreat_Jan2022.txt', MF_rates)
Z_FM_adj = inv_trans_surf('FM_surface_midpoints_MAP_specTreat_Jan2022.txt', FM_rates)

## Mean surface
Z_MF_mean_adj = inv_trans_surf('MF_surface_midpoints_mean_specTreat_Jan2022.txt', MF_rates)
Z_FM_mean_adj = inv_trans_surf('FM_surface_midpoints_mean_specTreat_Jan2022.txt', FM_rates)



#####
# 01/15/2022
# make previous spaghetti plots 
# maybe change the color scheme
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
# 01/24/2021 update: add title for better differentiation

# 02/17/2023 updates for manuscript----
# (1) larger text 
# (2) larger, more readable legend text/labels
# (3) add a "stacked" to showcase marginal source distribution
plot_source_freq <- function(Z, recipient = 'F',
                             age_breaks = seq(15,50,by=1), 
                             normalize = TRUE,
                             ylimits = c(0,0.1),
                             analysisType = 'flexible',
                             line_position = 'identity'){
  # group age first
  if(max(age_breaks) != 50){
    age_breaks = c(age_breaks, 50)
  }
  Z = Z %>% mutate(MAge_group = get_age_bin(MAge, breaks = age_breaks),
                   FAge_group = get_age_bin(FAge, breaks = age_breaks))
  
  if(recipient=='F'){
    # take the mean density within each recipient age group
    Z = Z %>% group_by(FAge_group, MAge) %>% 
      summarise(mean_dens = mean(density)) %>%
      select(rec_group = FAge_group, source_age = MAge, 
             density=mean_dens)
    
    xlab = 'Male source age'
    colorlab = 'Female age: '
    #tit = sprintf('M->F transmissions, %s types', analysisType)
    tit = sprintf('M-->F transmissions')
  }else{
    # take the mean density within each recipient age group
    Z = Z %>% group_by(MAge_group, FAge) %>% 
      summarise(mean_dens = mean(density)) %>%
      select(rec_group = MAge_group, source_age = FAge, 
             density=mean_dens)
    
    xlab = 'Female source age'
    colorlab = 'Male age: '
    #tit = sprintf('F->M transmissions, %s types', analysisType)
    tit = sprintf('F-->M transmissions')
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
      geom_line(aes(color = rec_group), 
                position = line_position) +
      labs(y='Transmission frequency', x=xlab, 
           color=colorlab, title = tit) +
      scale_y_continuous(limits = ylimits) +
      #ggthemes::scale_color_tableau('Classic Cyclic')+
      scale_color_manual(values = as.vector(pals::ocean.phase(length(age_breaks))))+
      theme_bw(base_size = 16) +
      theme(legend.position = 'bottom', 
            legend.title = element_text(size=14),
            legend.text = element_text(size=12),
            legend.key.width = unit(0.5, 'cm'),
            legend.key.height = unit(0.1,'cm')) +
      guides(color = guide_legend(nrow=(length(age_breaks)-2) %/% 6 + 1,
                                  byrow=TRUE))
  )
  
  # return the age grouped data
  Z
}

## plotting
## with different color palette and diff age bands
# 1. Fixed alloc

ylims = c(0,0.15)
ag_br = seq(15, 50, by = 3)
at = 'fixed'

pdf('source_freq_fixedAlloc_Jan2022_titled.pdf',
    height=5, width=7)

plot_source_freq(Z_MF_adj, recipient = 'F', 
                 age_breaks = ag_br, ylimits = ylims,
                 analysisType = at)
plot_source_freq(Z_FM_adj, recipient = 'M', 
                 age_breaks = ag_br, ylimits = ylims,
                 analysisType = at)
plot_source_freq(Z_MF_mean_adj, recipient = 'F',
                 age_breaks = ag_br, ylimits = ylims,
                 analysisType = at)
plot_source_freq(Z_FM_mean_adj, recipient = 'M',
                 age_breaks = ag_br, ylimits = ylims,
                 analysisType = at)

dev.off()

# 2. flexible labels
at = 'flexible'
ylims = c(0,0.15)
ag_br = seq(15, 50, by = 3)

# pdf('source_freq_specTreat_Jan2022_titled.pdf',
#     height=5, width=7)

pdf('source_freq_specTreat_Jan2022_titled_updated.pdf',
    height=5, width=7)

plot_source_freq(Z_MF_adj, recipient = 'F', 
                 age_breaks = ag_br, ylimits = ylims,
                 analysisType = at)
plot_source_freq(Z_FM_adj, recipient = 'M', 
                 age_breaks = ag_br, ylimits = ylims,
                 analysisType = at)
plot_source_freq(Z_MF_mean_adj, recipient = 'F',
                 age_breaks = ag_br, ylimits = ylims,
                 analysisType = at)
# plot_source_freq(Z_FM_mean_adj, recipient = 'M',
#                  age_breaks = ag_br, ylimits = ylims)

## a little bit more smoothing on the last plot
Z = plot_source_freq(Z_FM_mean_adj, recipient = 'M',
                     age_breaks = ag_br, ylimits = ylims,
                     analysisType = at)

# Z[Z$rec_group == '48-50' & Z$source_age == 41.5, 3] = 0.043
# Z[Z$rec_group == '48-50' & Z$source_age == 42.5, 3] = 0.045
# # re-normalize and plot
# Z = Z %>% group_by(rec_group) %>%
#   mutate(norm_density = density/sum(density)) %>% 
#   ungroup() %>% 
#   select(rec_group, source_age, density=norm_density)
# ggplot(Z, aes(x=source_age, y=density)) +
#   geom_line(aes(color = rec_group)) +
#   labs(y='Transmission frequency', x='Female source age', 
#        color='Male Age') +
#   scale_y_continuous(limits = ylims) +
#   #ggthemes::scale_color_tableau('Classic Cyclic')+
#   scale_color_manual(values = as.vector(pals::ocean.phase(length(ag_br))))+
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'bottom', 
#         legend.title = element_text(size=8),
#         legend.text = element_text(size=6),
#         legend.key.width = unit(0.3, 'cm'),
#         legend.key.height = unit(0.1,'cm')) +
#   guides(color = guide_legend(nrow=(length(ag_br)-1) %/% 12 + 1,
#                               byrow=TRUE))

dev.off()

## 02/17/2023: add stacked plot to show marginal source distributions

ylims = c(0, 0.1)

pdf('source_freq_specTreat_Jan2022_stacked.pdf',
    height=5, width=7)

plot_source_freq(Z_MF_adj, recipient = 'F', 
                 normalize = FALSE,
                 age_breaks = ag_br, ylimits = ylims,
                 analysisType = at,
                 line_position = 'stack')
plot_source_freq(Z_FM_adj, recipient = 'M', 
                 normalize = FALSE,
                 age_breaks = ag_br, ylimits = ylims,
                 analysisType = at,
                 line_position = 'stack')
plot_source_freq(Z_MF_mean_adj, recipient = 'F',
                 normalize = FALSE,
                 age_breaks = ag_br, ylimits = ylims,
                 analysisType = at,
                 line_position = 'stack')
# plot_source_freq(Z_FM_mean_adj, recipient = 'M',
#                  age_breaks = ag_br, ylimits = ylims)

## a little bit more smoothing on the last plot
Z = plot_source_freq(Z_FM_mean_adj, recipient = 'M',
                     normalize = FALSE,
                     age_breaks = ag_br, ylimits = ylims,
                     analysisType = at,
                     line_position = 'stack')

dev.off()


#####
# stuff on male/female transmitter freq. boxplot (over posterior samples)

# load estimates first
freq_MF = read_csv('MF_sc_dens_young_women.csv')
glimpse(freq_MF)

freq_MF = freq_MF %>% 
  select(monte_carlo_id, age = sc_age, freq)

# load sampling rates by age group
sc_samp_rates = readRDS('source_samp_rates_age.rds')
M_SC_samps = sc_samp_rates$M_SC_samps %>% 
  mutate(age = age_group + 0.5) %>%
  select(monte_carlo_id, age, prob)

# join and compute
adj_freq_MF = left_join(freq_MF, M_SC_samps)
glimpse(adj_freq_MF)

adj_freq_MF = adj_freq_MF %>% 
  mutate(adj_freq = freq/prob)

### a unified function to do all above wrangling and plotting
adj_freq <- function(freq_file, sc_gen = 'M', 
                     #rec_age = 'young', 
                     ylim = c(0, 0.012),
                     adjust = TRUE, normalize=TRUE,
                     caption = TRUE){
  
  # get raw freqs
  freq = read_csv(freq_file) %>% 
    mutate(age = sc_age - 0.5) %>%
    select(monte_carlo_id, age, freq, rec_age_lb, rec_age_ub)
  
  # get recipient age ranges
  rec_lb = freq$rec_age_lb[1]
  rec_ub = freq$rec_age_ub[1]
  age_range = paste0(rec_lb,'-',rec_ub)
  
  freq = freq %>% select(monte_carlo_id, age, freq)
  
  # get rates by source gender
  rates_ls = readRDS('source_samp_rates_age.rds')
  if(sc_gen == 'M'){
    rates = rates_ls$M_SC_samps
  }else{
    rates = rates_ls$F_SC_samps
  }
  # set to age band middle point
  rates = rates %>% 
    select(monte_carlo_id, age = age_group, prob)
  
  # adjust by dividing the sampling rates
  if(adjust){
    adj_freq = left_join(freq, rates) %>% 
      mutate(adj_freq = freq/prob)
  }else{
    adj_freq = freq %>% mutate(adj_freq = freq)
  }
  
  
  # normalize so freqs sum to 1
  if(normalize){
    adj_freq = adj_freq %>% group_by(age) %>%
      mutate(norm_freq = adj_freq/sum(adj_freq)) %>%
      ungroup() %>%
      select(monte_carlo_id, age, adj_freq = norm_freq)
  }else{
    adj_freq = adj_freq %>% select(monte_carlo_id, age, adj_freq)
  }
  
  # plotting
  if(sc_gen == 'M'){
    tit = paste('For female recipients aged', age_range)
    xx = paste('male transmitter age')
  }else{
    tit = paste('For male recipients aged', age_range)
    xx = paste('female transmitter age')
  }
  
  # caption for point type allocation
  if(caption){
    fnames = stringr::str_split(freq_file, '_')[[1]]
    if(stringr::str_starts(fnames[length(fnames)], 'fixed')){
      capt = 'fixed point types'
    }else{
      capt = 'ambiguous point types'
    }
  }else{
    capt = NULL
  }

  
  ## pad a row with age=50 to allow ticker at 50
  adj_freq = rbind(adj_freq, c(NA, 50, NA))
  
  print(
  ggplot(adj_freq, aes(x=as.character(age), y=adj_freq)) + 
    geom_boxplot(outlier.shape = NA, coef = 0, 
                 position = position_nudge(x = 0.5)#,
                 #color = "#00BFC4"
                 ) +
    scale_y_continuous(limits = ylim)+
    scale_x_discrete(breaks = as.character(seq(15,50,by=5))) +
    labs(x = xx, y='transmission frequency',
         title = tit, caption = capt) +
    theme_bw(base_size = 14) +
    theme(legend.position = 'none')
  )
  
  adj_freq
}

## plotting

pdf('adj_freq_plots_by_source_age.pdf',
    height = 4, width = 8)

# M->F
adj_freq_MF = adj_freq(freq_file = 'MF_sc_dens_young_women3.csv', 
                       ylim = c(0.0025, 0.006),
                       adjust = TRUE, normalize = TRUE)

adj_freq_MF2 = adj_freq(freq_file = 'MF_sc_dens_young_women_fixedAlloc3.csv',
                        ylim = c(0.0025, 0.006),
                        adjust = TRUE, normalize = TRUE)


adj_freq_MF = adj_freq(freq_file = 'MF_sc_dens_young_women2426.csv', 
                       ylim = c(0.0025, 0.006),
                       adjust = TRUE, normalize = TRUE)

adj_freq_MF2 = adj_freq(freq_file = 'MF_sc_dens_young_women2426_fixedAlloc.csv',
                        ylim = c(0.0025, 0.006),
                        adjust = TRUE, normalize = TRUE)

# F->M
adj_freq_FM = adj_freq(freq_file = 'FM_sc_dens_middle_men30.csv', 
                       sc_gen = 'F',
                       ylim = c(0.0025, 0.006),
                       adjust = TRUE, normalize = TRUE)

adj_freq_MF2 = adj_freq(freq_file = 'FM_sc_dens_middle_men30_fixedAlloc.csv',
                        sc_gen = 'F',
                        ylim = c(0.0025, 0.006),
                        adjust = FALSE, normalize = FALSE)

dev.off()


#####
# 01/24/2022 update
# add new plots to show boxplots of source age distribution 
# by age bands of recipients

## another helpful function
## get the middle point of the age bin
get_age_bin_midpoint <- function(x, breaks = seq(15,50,by=1)){
  num_bin = length(breaks) - 1
  res = character(length(x))
  for(i in 1:num_bin){
    lb = breaks[i]
    ub = breaks[i+1]
    res[x>lb & x<ub] = mean(c(lb,ub))
  }
  res
}

# do simple thing using the long format surfaces used for spaghetti plots
## 1. for each recipient age (middle point), sample `nsamp` source ages using the frequency/density
## 2. group by recipient age, draw box plots out of the samples
## (should allow broader age bands)
source_age_boxplot <- function(Z, nsamp = 1000,
                               recipient = 'F',
                               age_breaks = seq(15,50,by=1),
                               analysisType = 'flexible'){
  # group age first
  if(max(age_breaks) != 50){
    age_breaks = c(age_breaks, 50)
  }
  Z = Z %>% mutate(MAge_group = get_age_bin(MAge, breaks = age_breaks),
                   FAge_group = get_age_bin(FAge, breaks = age_breaks))
  
  # get source age density by recipient age group
  if(recipient=='F'){
    # take the mean density within each recipient age group
    Z = Z %>% group_by(FAge_group, MAge) %>% 
      summarise(mean_dens = mean(density)) %>%
      select(rec_group = FAge_group, source_age = MAge, 
             density=mean_dens)
    
    xlab = 'Female recipient age'
    ylab = 'Male source age'
    #colorlab = 'Female Age'
    tit = sprintf('M->F transmissions, %s types', analysisType)
  }else{
    # take the mean density within each recipient age group
    Z = Z %>% group_by(MAge_group, FAge) %>% 
      summarise(mean_dens = mean(density)) %>%
      select(rec_group = MAge_group, source_age = FAge, 
             density=mean_dens)
    
    xlab = 'Male recipient age'
    ylab = 'Female source age'
    #colorlab = 'Male Age'
    tit = sprintf('F->M transmissions, %s types', analysisType)
  }
  
  # sample source ages within each rec group
  rec_groups = unique(Z$rec_group)
  sc_ages = sort(unique(Z$source_age))
  All_samps = NULL
  for(rg in rec_groups){
    age_dens = Z %>% filter(rec_group == rg) %>%
      arrange(source_age) %>% select(density) %>%
      pull()
    age_samps = sample(sc_ages, nsamp, replace = TRUE, prob = age_dens)
    All_samps = c(All_samps, age_samps)
  }
  plot_dat = data.frame(rec_group = rep(rec_groups, each=nsamp),
                        source_age = All_samps)
  
  ## add dummy role
  plot_dat = rbind(plot_dat, data.frame(rec_group = c('50-51'), 
                                        source_age = NA))
  
  # produce boxplot
  print(
    ggplot(plot_dat, aes(x=rec_group, y=source_age)) +
      geom_boxplot(outlier.size = 0.3,
                   position = position_nudge(x = 0.5)) + 
      labs(y=ylab, x=xlab, title = tit) +  
      scale_x_discrete(labels = as.character(age_breaks)) +
      scale_y_continuous(limits = c(15,50)) +
      theme_bw(base_size = 14)
  )
  
  # don't return data
}

# try it
## fixed types
pdf('adj_source_age_boxplots_fixed.pdf',
    height = 4, width = 8)
at = 'fixed'
source_age_boxplot(Z_MF_adj, recipient = 'F',
                   nsamp = 5000,
                   analysisType = at)
source_age_boxplot(Z_FM_adj, recipient = 'M', 
                   nsamp = 5000,
                   analysisType = at)
source_age_boxplot(Z_MF_mean_adj, recipient = 'F',
                   nsamp = 5000,
                   analysisType = at)
source_age_boxplot(Z_FM_mean_adj, recipient = 'M',
                   nsamp = 5000,
                   analysisType = at)
dev.off()

## Flexible types
pdf('adj_source_age_boxplots_flexible.pdf',
    height = 4, width = 8)
source_age_boxplot(Z_MF_adj, recipient = 'F',
                   nsamp = 5000,
                   analysisType = at)
source_age_boxplot(Z_FM_adj, recipient = 'M', 
                   nsamp = 5000,
                   analysisType = at)
source_age_boxplot(Z_MF_mean_adj, recipient = 'F',
                   nsamp = 5000,
                   analysisType = at)
source_age_boxplot(Z_FM_mean_adj, recipient = 'M',
                   nsamp = 5000,
                   analysisType = at)
dev.off()

