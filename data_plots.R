# 04/22/2022
# make plots for manuscript

# 06/08/2022 update
# combine points colored with post.probs WITH 2d density with HDI contour lines (from "Density_2d_CIs.R")

library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(wesanderson)

setwd('~/Documents/Research/HIV_transmission_flow/')
#dat = read_csv('Rakai_data_2022_with_type_freqs.csv') # data capped at L>=0.2

dat = read_csv('Rakai_data_Jan2022.csv') # raw data uncapped...
dat = dat %>% filter(!is.na(POSTERIOR_SCORE_LINKED),
                     !is.na(POSTERIOR_SCORE_MF),
                     !is.na(POSTERIOR_SCORE_FM))

# 1. Data demo plot for Section 2 ----------------

## 05/02/2022 update: color points differently
## darker color in Xi's paper
## lighter color: points also included in my paper
dat = dat %>% 
  mutate(paper = if_else(POSTERIOR_SCORE_LINKED > 0.6 & 
                           (POSTERIOR_SCORE_MF > 0.66 |POSTERIOR_SCORE_FM > 0.66),
                         'Xi', 'Bu'))
## a. age pairs scatter plot + histogram------
scatter1 = ggplot(dat, aes(x=MALE_AGE_AT_MID,
                           y= FEMALE_AGE_AT_MID,
                           color = paper)) +
  geom_point(size = 1.8) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_manual(values = c('grey70','grey30')) +
  labs(x='male age (mid-observation)', y = 'female age (mid-observation)') +
  geom_text(x=15, y=47, label = sprintf('previous: N=%s\nnow: N=%s', 
                                       sum(dat$paper=='Xi'),
                                       nrow(dat)),
            hjust = 0, size = 5)+
  theme_bw(base_size = 14) +
  theme(legend.position = 'none')

(hist_male = ggplot(dat, aes(x=MALE_AGE_AT_MID)) +
  geom_histogram(color = 'gray20', fill = 'gray90', bins = 25) +
  theme_void())

(hist_female = ggplot(dat, aes(x=FEMALE_AGE_AT_MID)) +
    geom_histogram(color = 'gray20', fill = 'gray90', bins = 25) +
    theme_void() + 
    coord_flip())


(age_scatter = hist_male + plot_spacer() +
  scatter1 + hist_female +
  plot_layout(ncol = 2, 
              nrow = 2, 
              widths = c(4, 1),
              heights = c(1, 4)))

## b. plots for linkage score distribution and direct score distribution
## version 1: histograms-------
## 05/02/2022 update: direction score >0.66 or <0.33 ----> MF or FM
##                    otherwise ----> nothing
# dat = dat %>% 
#   mutate(direction_type = if_else(POSTERIOR_SCORE_MF >= 0.5, 
#                                   'likely M->F',
#                                   'likely F->M'))
dat = dat %>% 
  mutate(direction_type = case_when(POSTERIOR_SCORE_MF >= 0.66 ~ 'likely M->F',
                                    POSTERIOR_SCORE_MF <= 0.33 ~ 'likely F->M',
                                    TRUE ~ 'uncertatin'))
hist_link = ggplot(dat, aes(x=POSTERIOR_SCORE_LINKED)) +
  geom_histogram(color = 'gray20', fill = 'gray90', bins = 25) +
  geom_vline(xintercept = 0.6, linetype = 2, size = 1.5)+
  scale_x_continuous(limits = c(0,1),
                     labels = scales::percent)+
  scale_y_continuous(limits = c(0,50)) +
  labs(x='linkage score')+
  theme_bw(base_size = 14)

hist_direction = ggplot(dat, aes(x=POSTERIOR_SCORE_MF)) +
  geom_histogram(color = 'gray20', bins = 25, aes(fill=direction_type)) +
  scale_x_continuous(limits = c(0,1),
                     labels = scales::percent)+
  scale_y_continuous(limits = c(0,40)) +
  scale_fill_manual(values = c(wes_palette("Moonrise3")[2:1],'grey80')) +
  labs(x='M->F direction score',fill='')+
  theme_bw(base_size = 14) +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=10),
        legend.key.size = unit(0.4, 'cm'),
        legend.margin = margin(t = -5, r = 0, b = 0, l = 0, unit = "pt"))

ggarrange(age_scatter, 
          ggarrange(hist_link, hist_direction, nrow = 2, labels = c('B','C')),
          ncol = 2,
          labels = 'A')


## version 2: dot in increasing values------
## !!! remove those weird observations (where MF + FM scores != 1)
dat = dat %>% 
  filter(POSTERIOR_SCORE_MF + POSTERIOR_SCORE_FM == 1) %>% # filtering here!
  arrange(POSTERIOR_SCORE_LINKED) %>% 
  mutate(linkage_order = seq_along(POSTERIOR_SCORE_LINKED),
         direction_score = if_else(POSTERIOR_SCORE_MF >= 0.5, 
                                   POSTERIOR_SCORE_MF,
                                   POSTERIOR_SCORE_FM))
dots_link = ggplot(dat) +
  geom_point(color = "gray10", size = 1.2, aes(x=linkage_order, y=POSTERIOR_SCORE_LINKED)) +
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent) +
  #scale_x_continuous(breaks = c(0:5)*100) +
  scale_x_continuous(breaks = NULL, labels = NULL)+
  labs(y='linkage score',x='')+
  theme_bw(base_size = 14)

dots_direction = ggplot(dat) +
  geom_point(size = 1.5, 
             aes(x=linkage_order, 
                 #y=direction_score,
                 y = POSTERIOR_SCORE_MF,
                 color = direction_type)) +
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent) +
  #scale_x_continuous(breaks = c(0:5)*100) +
  scale_x_continuous(breaks = NULL, labels = NULL)+
  scale_color_manual(values = wes_palette("Moonrise3")[2:1]) +
  labs(y='direction score',x='pairs',color='')+
  theme_bw(base_size = 14)+
  theme(legend.position = 'bottom',
        legend.text = element_text(size=10),
        legend.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"))

ggarrange(age_scatter, 
          ggarrange(dots_link+rremove("xlab"), 
                    dots_direction, 
                    nrow = 2, 
                    heights = c(3,4),
                    labels = c('B','C')),
          ncol = 2,
          labels = 'A')


# 2. section 5: data colored with scores --------
dat2 = read_csv('Rakai_data_2022_with_type_freqs.csv') # using processed data in analysis

dat2Xi = dat2 %>% 
  filter(POSTERIOR_SCORE_LINKED >= 0.6) %>%
  mutate(direction = if_else(POSTERIOR_SCORE_MF >= 0.5, 
                             'M->F', 'F->M'))

## a. age pairs with pre-classification------
preMF  = ggplot(dat2Xi %>% 
                  filter(direction == 'M->F'), 
                aes(x=MALE_AGE_AT_MID,
                    y= FEMALE_AGE_AT_MID)) +
  geom_point(color = "gray30", size = 1.8) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  labs(x='male age', y = 'female age',
       title = 'With pre-classification') +
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20))

preFM = ggplot(dat2Xi %>% 
                 filter(direction == 'F->M'), 
               aes(x=MALE_AGE_AT_MID,
                   y= FEMALE_AGE_AT_MID)) +
  geom_point(color = "gray30", size = 1.8) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  labs(x='male age', y = 'female age') +
  theme_bw(base_size = 14)

## b. all age pairs with coloring ------
## version 1: color by posterior probs------
allMF = ggplot(dat2,
               aes(x=MALE_AGE_AT_MID,
                   y= FEMALE_AGE_AT_MID)) +
  geom_point(size = 1.8, aes(color = freq_MF)) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_distiller(type = "seq",
                        direction = 1,
                        palette = "Greys",
                        limits = c(0,1.0),
                        breaks = c(0,0.25,0.5,0.75,1.0),
                        labels = scales::percent(c(0,0.25,0.5,0.75,1.0)))+
  labs(x='male age', y = 'female age', 
       color='posterior\nM->F\nprobability',
       title = 'Including all data') +
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20))

allFM = ggplot(dat2,
               aes(x=MALE_AGE_AT_MID,
                   y= FEMALE_AGE_AT_MID)) +
  geom_point(size = 1.8, aes(color = freq_FM)) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_distiller(type = "seq",
                        direction = 1,
                        palette = "Greys",
                        limits = c(0,1.0),
                        breaks = c(0,0.25,0.5,0.75,1.0),
                        labels = scales::percent(c(0,0.25,0.5,0.75,1.0)))+
  labs(x='male age', y = 'female age', 
       color='posterior\nF->M\nprobability') +
  theme_bw(base_size = 14)

ggarrange(ggarrange(preMF+rremove("xlab")+rremove("x.text")+rremove('x.ticks'),
                    allMF+rremove("xylab")+rremove("axis.text")+rremove('ticks'), 
                    ncol=2,
                    widths = c(3.3,3.8)),
          ggarrange(preFM, 
                    allFM+rremove("ylab")+rremove("y.text")+rremove('y.ticks'), 
                    ncol=2,
                    widths = c(3.3,3.8)),
          nrow = 2,
          heights = c(3,3.2),
          labels = c('A','B'))

## version 2: color by direction scores------

allMFdirect = ggplot(dat2,
               aes(x=MALE_AGE_AT_MID,
                   y= FEMALE_AGE_AT_MID)) +
  geom_point(size = 1.8, aes(color = POSTERIOR_SCORE_MF)) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_distiller(type = "seq",
                        direction = 1,
                        palette = "Greys",
                        limits = c(0,1.0),
                        breaks = c(0,0.25,0.5,0.75,1.0),
                        labels = scales::percent(c(0,0.25,0.5,0.75,1.0)))+
  labs(x='male age', y = 'female age', 
       color='M->F\ndirection\nscore',
       title = 'Including all data') +
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20))

allFMdirect = ggplot(dat2,
               aes(x=MALE_AGE_AT_MID,
                   y= FEMALE_AGE_AT_MID)) +
  geom_point(size = 1.8, aes(color = POSTERIOR_SCORE_FM)) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_distiller(type = "seq",
                        direction = 1,
                        palette = "Greys",
                        limits = c(0,1.0),
                        breaks = c(0,0.25,0.5,0.75,1.0),
                        labels = scales::percent(c(0,0.25,0.5,0.75,1.0)))+
  labs(x='male age', y = 'female age', 
       color='F->M\ndirection\nscore') +
  theme_bw(base_size = 14)


ggarrange(ggarrange(preMF+rremove("xlab")+rremove("x.text")+rremove('x.ticks'),
                    allMFdirect+rremove("xylab")+rremove("axis.text")+rremove('ticks'), 
                    ncol=2,
                    widths = c(3.3,3.75)),
          ggarrange(preFM, 
                    allFMdirect+rremove("ylab")+rremove("y.text")+rremove('y.ticks'), 
                    ncol=2, 
                    widths = c(3.3,3.75)),
          nrow = 2,
          heights = c(3,3.2),
          labels = c('A','B'))


## 3. presentation illustration------
## contrasting fixed threshold approach and flexible approach
dat2 = read_csv('Rakai_data_2022_with_type_freqs.csv') # using processed data in analysis

dat2Xi = dat2 %>% 
  filter(POSTERIOR_SCORE_LINKED >= 0.6) %>%
  mutate(direction = case_when(POSTERIOR_SCORE_MF >= 0.66 ~ 'M->F',
                               POSTERIOR_SCORE_MF <= 0.33 ~ 'F->M',
                               TRUE ~ 'none'))

## 06/15: get data for fixed analysis as well
dat2fixed = dat2 %>%
  mutate(direction = case_when((POSTERIOR_SCORE_MF >= 0.5) & (POSTERIOR_SCORE_LINKED >= 0.6)  ~ 'M->F',
                               (POSTERIOR_SCORE_MF < 0.5) & (POSTERIOR_SCORE_LINKED >= 0.6) ~ 'F->M',
                               TRUE ~ 'none'))

## (1) Xi's model
## points
preMF  = ggplot(dat2Xi %>% 
                  filter(direction == 'M->F'), 
                aes(x=MALE_AGE_AT_MID,
                    y= FEMALE_AGE_AT_MID)) +
  geom_point(color = "gray30", size = 1.8) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  geom_text(x=20, y=48, label = sprintf('N=%s', sum(dat2Xi$direction=='M->F')),
            size = 12)+
  labs(x='male age', y = 'female age') +
  theme_bw(base_size = 14)


## discretized tiles
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
dat2XiMF = dat2Xi %>% 
  filter(direction == 'M->F') %>%
  mutate(male_age_group = get_age_bin(MALE_AGE_AT_MID),
         female_age_group = get_age_bin(FEMALE_AGE_AT_MID))

age_group_counts = dat2XiMF %>% 
  group_by(male_age_group, female_age_group) %>%
  summarise(num_points = n())

preMFtiles = ggplot(age_group_counts,
                    aes(x=male_age_group,
                        y=female_age_group,
                        fill = num_points)) +
  geom_tile() +
  scale_fill_gradient(low = alpha(wes_palette("Moonrise3")[1], 0.4), 
                      high = wes_palette("Darjeeling2")[2]) +
  labs(x='male age', y = 'female age') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



## (2) My model
## points
allMF = ggplot(dat2,
               aes(x=MALE_AGE_AT_MID,
                   y= FEMALE_AGE_AT_MID)) +
  geom_point(size = 1.8, aes(color = freq_MF)) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_gradient(low = 'grey85', high='black')+
  # scale_color_distiller(type = "seq",
  #                       direction = 1,
  #                       palette = "Greys",
  #                       limits = c(0,1.0),
  #                       breaks = c(0,0.25,0.5,0.75,1.0),
  #                       labels = scales::percent(c(0,0.25,0.5,0.75,1.0)))+
  geom_text(x=20, y=48, label = sprintf('N=%s', nrow(dat2)), size = 12)+
  labs(x='male age', y = 'female age', 
       color='posterior\nM->F\nprobability') +
  theme_bw(base_size = 14)+
  theme(legend.position = 'none')

## 2d density
allMFdens = ggplot(dat2 %>% filter(freq_MF > 0.6),
                   aes(x=MALE_AGE_AT_MID,
                       y= FEMALE_AGE_AT_MID)) +
  stat_density_2d(geom = "polygon", contour = TRUE, 
                  aes(fill = after_stat(level)), 
                  colour = "gray60",bins = 10) +
  # scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_fill_gradient(low = 'white', 
                      high = wes_palette("Darjeeling2")[2]) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  #geom_text(x=18, y=48, label = sprintf('N=%s', nrow(dat2)))+
  labs(x='male age', y = 'female age') +
  theme_bw(base_size = 14)+
  theme(legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



## 4. spatial points colored with posterior type probs-------
## the figure used in Section 5 of the paper

## (a) Xi paper classification
preMFcolored  = ggplot(dat2Xi %>% 
                  filter(direction == 'M->F'), 
                aes(x=MALE_AGE_AT_MID,
                    y= FEMALE_AGE_AT_MID)) +
  geom_point(color = wes_palette("Moonrise3")[1], size = 1.8) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  geom_text(x=18, y=48, label = sprintf('N=%s', sum(dat2Xi$direction=='M->F')), size = 6)+
  labs(x='', y = 'female age', title = 'M->F transmission') +
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 22))

preFMcolored  = ggplot(dat2Xi %>% 
                         filter(direction == 'F->M'), 
                       aes(x=MALE_AGE_AT_MID,
                           y= FEMALE_AGE_AT_MID)) +
  geom_point(color = wes_palette("Moonrise3")[2], size = 1.8) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  geom_text(x=18, y=48, label = sprintf('N=%s', sum(dat2Xi$direction=='F->M')), size = 6)+
  labs(x='', y = '', title = 'F->M transmission') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5,
                                   size = 22))

pre0colored = ggplot(dat2Xi %>% 
                       filter(direction == 'none'), 
                     aes(x=MALE_AGE_AT_MID,
                         y= FEMALE_AGE_AT_MID)) +
  geom_point(color = 'grey50', size = 1.8) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  geom_text(x=18, y=48, label = sprintf('N=%s', sum(dat2Xi$direction=='none')), size = 6)+
  labs(x='', y = '', title = 'no transmission link') +
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 22))


## (b) my paper flexible classification
# FM Pink: wes_palette("GrandBudapest1")[2]
# MF Blue: wes_palette("Darjeeling2")[2]

allMFcolored = ggplot(dat2,
               aes(x=MALE_AGE_AT_MID,
                   y= FEMALE_AGE_AT_MID)) +
  geom_point(size = 1.8, aes(color = freq_MF)) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_gradient(low = 'white', high=wes_palette("Darjeeling2")[2])+
  geom_text(x=18, y=48, label = sprintf('N=%s', nrow(dat2)), size = 6)+
  labs(x='male age', y = 'female age', 
       color='posterior\nM->F\nprobability') +
  theme_bw(base_size = 14)+
  theme(legend.position = 'none')


allFMcolored = ggplot(dat2,
               aes(x=MALE_AGE_AT_MID,
                   y= FEMALE_AGE_AT_MID)) +
  geom_point(size = 1.8, aes(color = freq_FM)) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_gradient(low = 'white', high=wes_palette("GrandBudapest1")[2])+
  geom_text(x=18, y=48, label = sprintf('N=%s', nrow(dat2)), size = 6)+
  labs(x='male age', y = '', 
       color='posterior\nF->M\nprobability') +
  theme_bw(base_size = 14)+
  theme(legend.position = 'none')

all0colored = ggplot(dat2,
                     aes(x=MALE_AGE_AT_MID,
                         y= FEMALE_AGE_AT_MID)) +
  geom_point(size = 1.8, aes(color = freq_0)) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_gradient(low = 'white', high='black')+
  geom_text(x=18, y=48, label = sprintf('N=%s', nrow(dat2)), size = 6)+
  labs(x='male age', y = '', 
       color='') +
  theme_bw(base_size = 14)+
  theme(legend.position = 'none')



ggarrange(ggarrange(preMFcolored,
                    preFMcolored,
                    pre0colored,
                    ncol=3),
          ggarrange(allMFcolored, 
                    allFMcolored,
                    all0colored,
                    ncol=3),
          nrow = 2,
          #heights = c(3,3.2),
          labels = c('A','B'))


# 6. colored points with spatial density contour lines as well------
## (June 11 updated version)

# 6.1 try the full model ones (only MF and FM directions) for proof of concept
# below code is copied over from "Density_2d_CIs.R" file

## the level function shared by Melodie M. 
getLevel <- function(x,y,z,prob=0.95) {
  dx <- diff(unique(x)[1:2])
  dy <- diff(unique(y)[1:2])
  sz <- sort(z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

## load addtional packages as needed
library(viridis)
library(wesanderson)

### try out additional color scales (for the points and the contours)
library(ggnewscale)

## look at those probs levels
probs = c(0.5, 0.8, 0.9)

# 1. MF surface first...
## the data for MF densities from full model
MFdens = read_csv('MF_density_MAP.csv')

pal = wes_palette("Moonrise2", 3, type='continuous')
MFlevels = MFdens %>% 
  summarise(level = getLevel(x,y,density, prob = probs)) %>%
  mutate(prob = probs)

MFlabels = MFdens %>% 
  full_join(MFlevels,by = character()) %>%
  mutate(diffs = abs(density - level)) %>%
  group_by(prob) %>%
  filter(diffs == min(diffs)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(prob_label = sprintf('%.0f%%', prob * 100))

## try some different text label positions here
## to avoid text-point overlap
MFlabels = MFdens %>% 
  full_join(MFlevels,by = character()) %>%
  mutate(diffs = abs(density - level)) %>%
  group_by(prob) %>%
  arrange(diffs) %>%
  slice(40) %>% # probably need to toggle this line and the next to find a good one
  filter(x == max(x)) %>%
  ungroup() %>%
  mutate(prob_label = sprintf('%.0f%%', prob * 100))

# reproduce "allMFcolored" but with contour lines of HDIs
allMFcolored_contours = 
  ggplot() +
  geom_point(data = dat2, 
             aes(x=MALE_AGE_AT_MID,
                 y= FEMALE_AGE_AT_MID,
                 color = freq_MF)) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_gradient(low = 'white', high=wes_palette("Darjeeling2")[2])+
  geom_text(data = dat2, 
            aes(x=MALE_AGE_AT_MID,
                y= FEMALE_AGE_AT_MID),
            x=18, y=48, label = sprintf('N=%s', nrow(dat2)), size = 6)+
  new_scale("color") +
  geom_contour(data = MFdens,
               aes(x=x,y=y, z=density, col = ..level..),
               breaks = MFlevels$level) +
  geom_text(data = MFlabels,
            aes(x=x, y=y, label = prob_label, col = level),
            size = 4) +
  scale_color_gradientn(colors = pal)+
  labs(x='male age', y = 'female age', 
       color='posterior\nM->F\nprobability') +
  theme_bw(base_size = 14)+
  theme(legend.position = 'none')


# 2. FM surface
FMdens = read_csv('FM_density_MAP.csv')

pal <- wes_palette("Zissou1", 3, type = "continuous")

FMlevels = FMdens %>% 
  summarise(level = getLevel(x,y,density, prob = probs)) %>%
  mutate(prob = probs)

FMlabels = FMdens %>% 
  full_join(FMlevels,by = character()) %>%
  mutate(diffs = abs(density - level)) %>%
  group_by(prob) %>%
  filter(diffs == min(diffs)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(prob_label = sprintf('%.0f%%', prob * 100))

## also try another way of labelling... 
## this one does NOT work that well!!
FMlabels = FMdens %>% 
  full_join(FMlevels,by = character()) %>%
  mutate(diffs = abs(density - level)) %>%
  group_by(prob) %>%
  arrange(diffs) %>%
  slice(30) %>% # probably need to toggle this line and the next to find a good one
  filter(x == max(x)) %>%
  ungroup() %>%
  mutate(prob_label = sprintf('%.0f%%', prob * 100))


## reproduce "allFMcolored" but with contour lines of HDIs
(
allFMcolored_contours = 
  ggplot() +
  geom_point(data = dat2, aes(color = freq_FM, 
                              x=MALE_AGE_AT_MID,
                              y= FEMALE_AGE_AT_MID), 
             size = 1.8) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  scale_color_gradient(low = 'white', high=wes_palette("GrandBudapest1")[2])+
  geom_text(data = dat2, aes(x=MALE_AGE_AT_MID,
                             y= FEMALE_AGE_AT_MID), 
            x=18, y=48, label = sprintf('N=%s', nrow(dat2)), size = 6)+
  new_scale("color") +
  geom_contour(data = FMdens,
               aes(x=x,y=y, z=density, col = ..level..),
               breaks = FMlevels$level) +
  geom_text(data = FMlabels,
            aes(x=x, y=y, label = prob_label, col = level),
            size = 4) +
  scale_color_gradientn(colors = pal)+
  labs(x='male age', y = 'female age', 
       color='posterior\nF->M\nprobability') +
  theme_bw(base_size = 14)+
  theme(legend.position = 'none')
)


## 6.2 then also the fixed analysis ones----

## (a) MF surface

## the data for MF densities from partial model
MFdens = read_csv('fixThres_MF_density_MAP.csv')

pal = wes_palette("Moonrise2", 3, type='continuous')
MFlevels = MFdens %>% 
  summarise(level = getLevel(x,y,density, prob = probs)) %>%
  mutate(prob = probs)

MFlabels = MFdens %>% 
  full_join(MFlevels,by = character()) %>%
  mutate(diffs = abs(density - level)) %>%
  group_by(prob) %>%
  filter(diffs == min(diffs)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(prob_label = sprintf('%.0f%%', prob * 100))

## try some different text label positions here
## to avoid text-point overlap
MFlabels = MFdens %>% 
  full_join(MFlevels,by = character()) %>%
  mutate(diffs = abs(density - level)) %>%
  group_by(prob) %>%
  arrange(diffs) %>%
  slice(40) %>% # probably need to toggle this line and the next to find a good one
  filter(x == max(x)) %>%
  ungroup() %>%
  mutate(prob_label = sprintf('%.0f%%', prob * 100))

(
preMFcolored  = 
  ggplot() +
  geom_point(data = dat2fixed %>% 
               filter(direction == 'M->F'), 
             aes(x=MALE_AGE_AT_MID,
                 y= FEMALE_AGE_AT_MID),
             color = wes_palette("Moonrise3")[1], size = 1.8) +
  scale_x_continuous(limits = c(15,50)) +
  scale_y_continuous(limits = c(15,50)) +
  geom_text(data = dat2fixed, 
            aes(x=MALE_AGE_AT_MID,
                y= FEMALE_AGE_AT_MID),x=18, y=48, 
            label = sprintf('N=%s', sum(dat2fixed$direction=='M->F')), size = 6)+
  geom_contour(data = MFdens,
               aes(x=x,y=y, z=density, col = ..level..),
               breaks = MFlevels$level) +
  geom_text(data = MFlabels,
            aes(x=x, y=y, label = prob_label, col = level),
            size = 4) +
  scale_color_gradientn(colors = pal)+
  labs(x='male age', y = 'female age', title = 'M->F transmissions') +
  theme_bw(base_size = 14)+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5,
                                          size = 22))
)


## (b) FM surface
FMdens = read_csv('fixThres_FM_density_MAP.csv')

pal <- wes_palette("Zissou1", 3, type = "continuous")

FMlevels = FMdens %>% 
  summarise(level = getLevel(x,y,density, prob = probs)) %>%
  mutate(prob = probs)

FMlabels = FMdens %>% 
  full_join(FMlevels,by = character()) %>%
  mutate(diffs = abs(density - level)) %>%
  group_by(prob) %>%
  filter(diffs == min(diffs)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(prob_label = sprintf('%.0f%%', prob * 100))

## also try another way of labelling... 
## this one does NOT work that well!!
FMlabels = FMdens %>% 
  full_join(FMlevels,by = character()) %>%
  mutate(diffs = abs(density - level)) %>%
  group_by(prob) %>%
  arrange(diffs) %>%
  slice(30) %>% # probably need to toggle this line and the next to find a good one
  filter(x == max(x)) %>%
  ungroup() %>%
  mutate(prob_label = sprintf('%.0f%%', prob * 100))


(
  preFMcolored  = 
    ggplot() +
    geom_point(data = dat2fixed %>% 
                 filter(direction == 'F->M'), 
               aes(x=MALE_AGE_AT_MID,
                   y= FEMALE_AGE_AT_MID),
               color = wes_palette("Moonrise3")[2], size = 1.8) +
    scale_x_continuous(limits = c(15,50)) +
    scale_y_continuous(limits = c(15,50)) +
    geom_text(data = dat2fixed, 
              aes(x=MALE_AGE_AT_MID,
                  y= FEMALE_AGE_AT_MID),x=18, y=48, 
              label = sprintf('N=%s', sum(dat2fixed$direction=='F->M')), size = 6)+
    geom_contour(data = FMdens,
                 aes(x=x,y=y, z=density, col = ..level..),
                 breaks = FMlevels$level) +
    geom_text(data = FMlabels,
              aes(x=x, y=y, label = prob_label, col = level),
              size = 4) +
    scale_color_gradientn(colors = pal)+
    labs(x='male age', y = 'female age', title = 'F->M transmissions') +
    theme_bw(base_size = 14)+
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5,
                                    size = 22))
)


## put the four panels (colored points + density contour lines) together
ggarrange(ggarrange(preMFcolored,
                    preFMcolored,
                    ncol=2),
          ggarrange(allMFcolored_contours, 
                    allFMcolored_contours,
                    ncol=2),
          nrow = 2,
          heights = c(3.2,3),
          labels = c('A','B'))
