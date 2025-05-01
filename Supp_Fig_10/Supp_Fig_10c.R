library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Load metadata paCoV cohort
Meta <- read.csv('./Meta_paCoV40_03112025_small.csv')

# colors 
col_patient_groups = c('aG1'= '#deebf7', 'aG2'='#9ecae1', 'aG3'='#4292c6', 
                       'aHC'='#addd8e', 'pG1'='#d8daeb', 'pG2'='#9e9ac8', 'pG3'='#54278f', 'pHC'='#66bd63')

# Order Subclusters
SC_order <- c('B_cells_SC0', 'B_cells_SC1', 'B_cells_SC2', 'B_cells_SC3',
              'B_cells_SC4', 'B_cells_SC5', 'B_cells_SC6', 'CD4_SC0', 'CD4_SC1',
              'CD4_SC2', 'CD4_SC3', 'CD4_SC4', 'CD4_SC5', 'CD8_SC0', 'CD8_SC1',
              'CD8_SC2', 'CD8_SC3', 'CD8_SC4', 'CD14_SC0', 'CD14_SC1', 'CD14_SC2',
              'CD14_SC3', 'CD16_mono_SC0', 'CD16_mono_SC1', 'CD16_mono_SC2',
              'CD16_mono_SC3', 'Eryth', 'HSCs', 'Mgk_SC0', 'Mgk_SC1', 'Mgk_SC2',
              'Mgk_SC3', 'NK_SC0', 'NK_SC1', 'PC_SC0', 'PC_SC1', 'PC_SC2', 'PC_SC3',
              'PC_SC4', 'PC_SC5', 'PC_SC6', 'cDC_SC0', 'cDC_SC1', 'cDC_SC2',
              'cDC_SC3', 'pDC_SC0', 'pDC_SC1', 'pDC_SC2')


BP <- Meta %>% 
  mutate(Patient_groups = factor(Patient_groups, levels = c("pHC","pG1","pG2",'pG3',
                                                            "aHC","aG1","aG2",'aG3'))) %>%
  mutate(ReCluster = factor(SCs)) %>%
  group_by( ReCluster, Patient_groups) %>%
  summarise(n = n()) %>% #, Set = first(Set)
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>%
  ggplot(aes(x = ReCluster, y = freq, fill = Patient_groups, group = Patient_groups)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=col_patient_groups) +
  scale_x_discrete(limits=SC_order) + 
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle = 90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") + #    ylab('% PBMC') + xlab('Age groups')
  
  ylab('Prop. of cells') + xlab('Individuals (n=40)')

print(BP)
