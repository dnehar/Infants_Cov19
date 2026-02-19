library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# colors
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# load meta data (PBMCs, n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')
head(Meta)
dim(Meta) 

clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups,2, FUN = list, simplify = T)

# subsets to plot
subset_to_be_plotted <- c( 'CD4_SC0','CD4_SC1','CD4_SC2','CD4_SC3')

# plot
plt_clinical <- Meta %>% 
  mutate(ReCluster = factor(annotated_SCs)) %>%
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  filter(ReCluster %in% subset_to_be_plotted) %>%
  group_by(Groups, Fig_ids, ReCluster) %>%
  summarise(n = n()) %>% #, Set = first(Set)
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>% 
  #filter(ReCluster %in% subset_to_be_plotted) %>% 
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.2) +
  theme_bw()  +  #THEME +
  #ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  #ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = F, vjust = 0.5) + 
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face='bold')) +
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  
  scale_fill_manual(values=col_pGroups) + #**
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =0),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) + #    ylab('% PBMC') + xlab('Age groups')
  ylab('% in CD4 T cells') + xlab('Clinical groups')

print(plt_clinical)

