library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)


# load meta data (PBMCs, n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')
head(Meta)
dim(Meta) 

# colors
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups,2, FUN = list, simplify = T)

BP <- Meta %>% 
  mutate(Patient_groups = factor(Patient_groups, levels = c("pHC", "G1", "G2", "G3"))) %>%
  mutate(ReCluster = factor(annotated_SCs)) %>%
  group_by( ReCluster, Patient_groups) %>%
  #filter(Groups %in% c("HO_M",'HO_F')) %>% 
  summarise(n = n()) %>% #, Set = first(Set)
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>%
  ggplot(aes(x = ReCluster, y = freq, fill = Patient_groups, group = Patient_groups)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=col_pGroups) + #***
  #scale_x_discrete(limits=c("pHC", "G1", "G2", "G3")) + #labels= labels
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle = 90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") + #    ylab('% PBMC') + xlab('Age groups')
  
  ylab('Prop. of cells') + xlab('Individuals (n=40)')

print(BP) 
