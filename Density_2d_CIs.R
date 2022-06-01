# 05/30/2022

# Use Melodie's code to compute 2d credible regions for inference results (the density surfaces)

# her function
getLevel <- function(x,y,z,prob=0.95) {
  dx <- diff(unique(x)[1:2])
  dy <- diff(unique(y)[1:2])
  sz <- sort(z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

## demo code ------
library(data.table)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(viridis)

set.seed(1001)
d <- data.table(expand.grid(x=seq(-4,4, 0.1),y = seq(-4,4, 0.1))) # your domain
d[, density := dmvnorm(c(x,y), mean = c(0,0), sigma = matrix(c(1, 3/5, 3/5, 1), nrow =2, byrow = T)), by = c('x', 'y')] # this is your spatial density

levels <- d[, {
  prob = c(0.5, 0.8, 0.9)
  level = getLevel(x, y, density, prob)
  list(prob = prob, level = level)
}]

d = as.data.table(full_join(d, levels, by = character()))

# find where to place the labels
d[, diff := abs(density - level)]
mindiff <- d[, list(diff = min(diff)), by = 'prob']
d1 <- merge(d, mindiff, by = c('prob', 'diff'))
d1 <- d1[, list(x = x[1], y = y[1]), by = c('prob', 'level')] # pick the first coordinates
d1[, prob_label := paste0(prob * 100, '%')]

ggplot(d,aes(x=x,y=y))+
  geom_raster(aes(fill = density)) +
  geom_contour(aes(z=density, col = ..level..), breaks = levels[, level]) +
  geom_text(data = d1, aes(label = prob_label, col = level), size = 4) +
  scale_color_viridis_c() +
  scale_fill_gradient(low = 'beige', high = 'firebrick') +
  guides(col="none")

## demo code above ##



## my code ------

library(wesanderson)

probs = c(0.5, 0.8, 0.9)

### 1. MF surface
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

p_MF = ggplot(MFdens,aes(x=x,y=y))+
  geom_raster(aes(fill = density)) +
  geom_contour(aes(z=density, col = ..level..), breaks = MFlevels$level) +
  geom_text(data = MFlabels, aes(label = prob_label, col = level), size = 4) +
  annotate(geom = 'text', x=21, y=47, 
           label = 'M->F transmission', size = 6)+
  #scale_color_viridis_c(option = 'C') +
  scale_color_gradientn(colors = pal)+
  scale_fill_gradient(low = 'white', 
                      high = wes_palette("Darjeeling2")[2]) +
  scale_x_continuous('male age', expand = c(0, 0)) + 
  scale_y_continuous('female age', expand = c(0, 0)) +
  guides(col="none", fill='none')+
  theme_bw(base_size = 14)

### 2. FM surface
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

p_FM = 
ggplot(FMdens,aes(x=x,y=y))+
  geom_raster(aes(fill = density)) +
  geom_contour(aes(z=density, col = ..level..), breaks = FMlevels$level) +
  geom_text(data = FMlabels, aes(label = prob_label, col = level), size = 4) +
  annotate(geom = 'text', x=21, y=47, 
           label = 'F->M transmission', size = 6)+
  #scale_color_viridis_c(option = 'E') +
  scale_color_gradientn(colors = pal)+
  scale_fill_gradient(low = 'white', 
                      high = wes_palette("GrandBudapest1")[2]) +
  scale_x_continuous('male age', expand = c(0, 0)) + 
  scale_y_continuous('female age', expand = c(0, 0)) +
  guides(col="none", fill='none')+
  theme_bw(base_size = 14)
