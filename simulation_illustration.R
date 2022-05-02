# 04/24/2022
# try producing an illutration for simulations setup

library(dplyr)
library(ggplot2)

setups = data.frame(prop = c(0.3, 0.6, 0.1,  0.6, 0.3, 0.1, 0.5, 0.5, 0.6, 0.4),
                    type = c('men ~25y.o.', 'men ~35y.o.', 'else', 'men ~25y.o.', 'men ~35y.o.', 'else',
                             'M->F', 'F->M', 'M->F', 'F->M'),
                    scenario = c(rep('discordant age',3), 
                                 rep('same age', 3),
                                 rep('MF 50-50', 2),
                                 rep('MF 60-40', 2)),
                    problem = c(rep('age',6), rep('direction',4)))
labels = data.frame(label = c('men ~35y.o.\n60%','men ~25y.o.\n30%',  'else 10%', 
                              'men ~35y.o.\n30%', 'men ~25y.o.\n60%', 'else 10%',
                              'M->F\n50%', 'F->M\n50%', 'M->F\n60%', 'F->M\n40%'),
                    scenario = c(rep('discordant age',3), 
                                 rep('same age', 3),
                                 rep('MF 50-50', 2),
                                 rep('MF 60-40', 2)),
                    y = c(0.3, 0.75, 0.95,
                          0.15, 0.6, 0.95, 
                          0.25, 0.75,
                          0.35, 0.8))


ggplot(setups, aes(x=scenario)) +
  geom_bar(stat='identity', aes(y=prop, fill=type), alpha = 0.9) +
  geom_label(data=labels,aes(label=label,x=scenario, y=y))+
  geom_vline(xintercept = 2.5, size = 1.2)+
  scale_fill_manual(breaks = c('men ~25y.o.', 'men ~35y.o.', 'F->M', 'M->F', 'else'),
                      values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF", "gray80")) +
  scale_x_discrete(breaks = c('same age','discordant age','MF 50-50','MF 60-40'),
                   limits = c('same age','discordant age','MF 50-50','MF 60-40'),
                   labels = c('Scenario1\nsame age','Scenario2\ndiscordant age',
                              'Scenario1\nMF 50-50','Scenario2\nMF 60-40'),
                   position = 'top')+
  # scale_y_continuous(labels = scales::percent,
  #                    breaks = c(0,0.3,0.5,0.6,0.9,1.0))+
  scale_y_continuous(NULL, breaks = NULL)+
  labs(x='', y='source proportions')+
  theme_minimal(base_size = 16) +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=15,face='bold'))
