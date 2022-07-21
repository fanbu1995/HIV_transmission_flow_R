# 05/02/2022 update
# simple summary statistics on real analysis

# 05/30/2022 update
# add contour plots for MF/FM surface estimates

# we have:
# 1. proportions of each type
# 2. spatial points colored with posterior type probs
# 3. source age distribution (marginal on male or female)
# 4. density of MF and FM surface, with contours for HDIs (in file "Density_2d_CIs.R")
# 5. combine the colored points (with post.probs) with contour lines of HDIs (in "data_plot.R")


library(tidyverse)
library(bayestestR)
library(xtable)
setwd('~/Documents/Research/HIV_transmission_flow/')

# 1. proportions of each type --------
probs = read_csv('real_data_probs_chains.csv')
fixed_probs = read_csv('fixed_real_data_probs_chains.csv')

burn = 1000

probs_summary = apply(probs[(burn+1):nrow(probs),], 2, 
                      function(v) c(mean(v, na.rm = TRUE),
                                    quantile(v, c(0.025, 0.975), na.rm = TRUE),
                                    sd(v, na.rm = TRUE)))

probs_summary = as.data.frame(t(probs_summary))
names(probs_summary) = c('avg', 'CI95_lb', 'CI95_ub', 'std')
probs_summary$prob = names(probs)

fixed_summary = data.frame(prob = names(fixed_probs),
                           avg = fixed_probs[1,] %>% as.numeric())

## 06/29/2022: output a summary table as well
comb_summary = cbind(probs_summary %>% select(prob, avg, CI95_lb, CI95_ub),
                     fixed_summary %>% select(fixed_avg = avg)) %>%
  mutate(model_N = round(avg * 526, digits = 0), 
         fixed_N = fixed_avg * 526) %>%
  mutate(model_avg_text = 
           sprintf('%.1f%% (%.1f%%, %.1f%%)', 
                   avg * 100, CI95_lb * 100, CI95_ub * 100),
         fixed_avg_text = sprintf('%.1f%%', fixed_avg * 100)) %>%
  select(prob, model_avg_text, model_N,
         fixed_avg_text, fixed_N)

print(xtable(comb_summary[c(2,3,1),]), 
      include.rownames = FALSE)


text_info = data.frame(x=names(probs),
                       y=0.8,
                       label = sprintf('fixed N=%.0f\n model N=%.0f',
                                       fixed_summary$avg * 526,
                                       probs_summary$avg * 526))


ggplot(probs_summary, aes(x=prob, y=avg)) +
  geom_bar(stat='identity', aes(color = prob, fill = prob))+
  geom_errorbar(aes(ymax = CI95_ub, ymin = CI95_lb), color='grey30', width = 0.5) +
  geom_point(data=fixed_summary, shape = 17, aes(x=prob, y=avg), size = 4)+
  geom_text(data=text_info, aes(x=x, y=y, label=label), hjust=0.5, size = 5)+
  scale_y_continuous(limits = c(0,1), labels = scales::percent) +
  scale_x_discrete(labels = c('none', 'M->F', 'F->M'))+
  scale_fill_manual(values = c('grey70',wes_palette("Moonrise3")[1:2]))+
  scale_color_manual(values = c('grey70',wes_palette("Moonrise3")[1:2]))+
  labs(x='type', y='posterior mean & 95% CI')+
  theme_bw(base_size = 14) +
  theme(legend.position = 'none')





# 2. spatial points colored with posterior type probs-------
# SEE "data_plots.R" file



# 3. source age distribution (marginal on male or female) ------
# 05/03 update: add age distribution from the fixed threshold analysis as well
males = read_csv('real_data_male_source_age_samples.csv')
females = read_csv('real_data_female_source_age_samples.csv')

malesFixed = read_csv('fixedThres_data_male_source_age_samples.csv')
femalesFixed = read_csv('fixedThres_data_female_source_age_samples.csv')

burn = 1000
thin = 20
maxIter = max(males$iter)

sel = seq(from=burn, to = maxIter, by=thin)

males = males %>% filter(iter %in% sel)
females = females %>% filter(iter %in% sel)

xlims = c(15,50)

## male source age
hdi50male = hdi(males$sourceAge, ci=0.5)
hdi50maleFixed = hdi(malesFixed$sourceAge, ci=0.5)

ggplot(males, aes(x=sourceAge)) +
  geom_density(aes(group=iter),
               color = alpha(wes_palette("Moonrise3")[1],0.2))+
  stat_density(bw=2,color = wes_palette("Darjeeling2")[2], 
               geom='line', position='identity', size = 1.5)+
  stat_density(data=malesFixed, aes(x=sourceAge), 
               bw=2, color='grey30', geom='line', linetype = 2,
               position='identity', size = 1.5) +
  geom_text(x=40, y = 0.06,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50male$CI_low, hdi50male$CI_high),
            size = 6, color = wes_palette("Darjeeling2")[2])+
  geom_text(x=40, y = 0.056,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50maleFixed$CI_low, hdi50maleFixed$CI_high),
            size = 6)+
  scale_x_continuous(limits = xlims)+
  labs(x='male source age', y='')+
  theme_bw(base_size = 14)


## female source age
hdi50female = hdi(females$sourceAge, ci=0.5)
hdi50femaleFixed = hdi(femalesFixed$sourceAge, ci=0.5)
ggplot(females, aes(x=sourceAge)) +
  geom_density(aes(group=iter), 
               color = alpha(wes_palette("Moonrise3")[2],0.2))+
  stat_density(bw=2.2,color = wes_palette("GrandBudapest1")[2], 
               geom='line', position='identity', size = 1.5)+
  stat_density(data=femalesFixed, aes(x=sourceAge), 
               bw=2, color='grey30', geom='line', linetype = 2,
               position='identity', size = 1.5) +
  geom_text(x=40, y = 0.06,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50female$CI_low, hdi50female$CI_high),
            size = 6, color = wes_palette("GrandBudapest1")[2])+
  geom_text(x=40, y = 0.056,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50femaleFixed$CI_low, hdi50femaleFixed$CI_high),
            size = 6)+
  labs(x='female source age', y='')+
  theme_bw(base_size = 14)

## 3.b. source/recipient age distribution for specific age group----

## e.g. I: young women (15-25)
source.males = read_csv('real_data_young_women_infection_male_source_age_samples.csv')
rec.males = read_csv('real_data_young_female_recipient_age_samples.csv')

source.males.fixed = read_csv('fixedThres_data_young_women_infection_male_source_age_samples.csv')
rec.males.fixed = read_csv('fixedThres_data_young_female_recipient_age_samples.csv')

source.males = source.males %>% filter(iter %in% sel)
rec.males = rec.males %>% filter(iter %in% sel)

source.males.fixed = source.males.fixed %>% filter(iter %in% sel)
rec.males.fixed = rec.males.fixed %>% filter(iter %in% sel)

ylims = c(0, 0.13)
xlims = c(15,50)

## male source age
hdi50male = hdi(source.males$sourceAge, ci=0.5)
hdi50male.fixed = hdi(source.males.fixed$sourceAge, ci=0.5)
ggplot(source.males, aes(x=sourceAge)) +
  geom_density(aes(group=iter), 
               color = alpha(wes_palette("Moonrise3")[1],0.2))+
  stat_density(bw=2,color = wes_palette("Darjeeling2")[2], 
               geom='line', position='identity', size = 1.5)+
  stat_density(data=source.males.fixed, aes(x=sourceAge), 
               bw=2, color='grey30', geom='line', linetype = 2,
               position='identity', size = 1.5) +
  geom_text(x=40, y = 0.10,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50male$CI_low, hdi50male$CI_high),
            size = 6, color = wes_palette("Darjeeling2")[2])+
  geom_text(x=40, y = 0.092,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50male.fixed$CI_low, hdi50male.fixed$CI_high),
            size = 6)+
  scale_y_continuous(limits = ylims)+
  scale_x_continuous(limits = xlims)+
  labs(x='male source age for young women (15-25)', y='')+
  theme_bw(base_size = 14)

## male recipient age
hdi50male = hdi(rec.males$sourceAge, ci=0.5)
hdi50male.fixed = hdi(rec.males.fixed$sourceAge, ci=0.5)
ggplot(rec.males, aes(x=sourceAge)) +
  geom_density(aes(group=iter), 
               color = alpha(wes_palette("Moonrise3")[2],0.2))+
  stat_density(bw=1.8,color = wes_palette("GrandBudapest1")[2], 
               geom='line', position='identity', size = 1.5)+
  stat_density(data=rec.males.fixed, aes(x=sourceAge), 
               bw=1.8, color='grey30', geom='line', linetype = 2,
               position='identity', size = 1.5) +
  geom_text(x=40, y = 0.10,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50male$CI_low, hdi50male$CI_high),
            size = 6, color = wes_palette("GrandBudapest1")[2])+
  geom_text(x=40, y = 0.092,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50male.fixed$CI_low, hdi50male.fixed$CI_high),
            size = 6)+
  scale_y_continuous(limits = ylims)+
  scale_x_continuous(limits = xlims)+
  labs(x='male recipient age from young women (15-25)', y='')+
  theme_bw(base_size = 14)
  
  
## e.g. II: young men (15-25)
source.females = read_csv('real_data_young_men_infection_female_source_age_samples.csv')
rec.females = read_csv('real_data_young_male_recipient_age_samples.csv')

source.females = source.females %>% filter(iter %in% sel)
rec.females = rec.females %>% filter(iter %in% sel)

source.females.fixed = read_csv('fixedThres_data_young_men_infection_female_source_age_samples.csv')
rec.females.fixed = read_csv('fixedThres_data_young_male_recipient_age_samples.csv')

source.females.fixed = source.females.fixed %>% filter(iter %in% sel)
rec.females.fixed = rec.females.fixed %>% filter(iter %in% sel)

ylims = c(0, 0.13)
xlims = c(15,50)

## female source age
hdi50female = hdi(source.females$sourceAge, ci=0.5)
hdi50femaleFixed = hdi(source.females.fixed$sourceAge, ci=0.5)
ggplot(source.females, aes(x=sourceAge)) +
  geom_density(aes(group=iter), 
               color = alpha(wes_palette("Moonrise3")[2],0.2))+
  stat_density(bw=2.2,color = wes_palette("GrandBudapest1")[2], 
               geom='line', position='identity', size = 1.5)+
  stat_density(data=source.females.fixed, aes(x=sourceAge), 
               bw=2.2, color='grey30', geom='line', linetype = 2,
               position='identity', size = 1.5) +
  geom_text(x=40, y = 0.10,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50female$CI_low, hdi50female$CI_high),
            size = 6, color = wes_palette("GrandBudapest1")[2])+
  geom_text(x=40, y = 0.092,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50femaleFixed$CI_low, hdi50femaleFixed$CI_high),
            size = 6)+
  scale_y_continuous(limits = ylims)+
  scale_x_continuous(limits = xlims)+
  labs(x='female source age for young men (15-25)', y='')+
  theme_bw(base_size = 14)

## female recipient age
hdi50female = hdi(rec.females$sourceAge, ci=0.5)
hdi50femaleFixed = hdi(rec.females.fixed$sourceAge, ci=0.5)
ggplot(rec.females, aes(x=sourceAge)) +
  geom_density(aes(group=iter), 
               color = alpha(wes_palette("Moonrise3")[1],0.2))+
  stat_density(bw=2,color = wes_palette("Darjeeling2")[2], 
               geom='line', position='identity', size = 1.5)+
  stat_density(data=rec.females.fixed, aes(x=sourceAge), 
               bw=2, color='grey30', geom='line', linetype = 2,
               position='identity', size = 1.5) +
  geom_text(x=40, y = 0.10,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50female$CI_low, hdi50female$CI_high),
            size = 6, color = wes_palette("Darjeeling2")[2])+
  geom_text(x=40, y = 0.092,
            label = sprintf('50%% HDI:[%.1f, %.1f]', 
                            hdi50femaleFixed$CI_low, hdi50femaleFixed$CI_high),
            size = 6)+
  scale_y_continuous(limits = ylims)+
  scale_x_continuous(limits = xlims)+
  labs(x='female recipient age from young men (15-25)', y='')+
  theme_bw(base_size = 14)
  