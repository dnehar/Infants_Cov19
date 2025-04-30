library(Seurat)
library(reshape2)
library(pheatmap)
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# clinical groups 
clinical_groups <- c('pHC',"pG1", 'pG2', 'pG3',
                     'aHC',"aG1", 'aG2', 'aG3')

my_comp <- list(c('pG2','aG1'), c('pG2','aG2'), c('pG2','aG3'), 
                c('pG3','aG1'), c('pG3','aG2'), c('pG3','aG3'))

# colors 
col_patient_groups = c('aG1'= '#deebf7', 'aG2'='#9ecae1', 'aG3'='#4292c6', 
                       'aHC'='#addd8e', 'pG1'='#d8daeb', 'pG2'='#9e9ac8', 'pG3'='#54278f', 'pHC'='#66bd63')

# Load metadata 
Meta <- read.csv('./Meta_paCoV40_03112025_small.csv')

######################### CD14 monocytes subsets #########################

subset_to_be_plotted <- c( 'CD14_SC0','CD14_SC1','CD14_SC2','CD14_SC3')

plt_CD14_mo <- Meta %>% 
  mutate(ReCluster = factor(SCs)) %>%
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
  ggpubr::stat_compare_means(comparisons = my_comp, method = "t.test", label = "p.signif") +
  #ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = F, vjust = 0.5) + 
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face='bold')) +
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  
  scale_fill_manual(values=col_patient_groups) + #***
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) + #    ylab('% PBMC') + xlab('Age groups')
  ylab('% in lineage') + xlab('Clinical groups')

print(plt_CD14_mo)

######################### CD4 T cell subsets ########################

subset_to_be_plotted <- c('CD4_SC0','CD4_SC1','CD4_SC2','CD4_SC3','CD4_SC4','CD4_SC5')

plt_CD4_T <- Meta %>% 
  mutate(ReCluster = factor(SCs)) %>%
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
  ggpubr::stat_compare_means(comparisons = my_comp, method = "t.test", label = "p.signif") +
  #ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = F, vjust = 0.5) + 
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face='bold')) +
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  
  scale_fill_manual(values=col_patient_groups) + #***
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) + #    ylab('% PBMC') + xlab('Age groups')
  ylab('% in CD4 T cells') + xlab('Clinical groups')

print(plt_CD4_T)

########################## B cell subsets #########################

subset_to_be_plotted <- c( 'B_cells_SC0','B_cells_SC1','B_cells_SC2','B_cells_SC3','B_cells_SC4','B_cells_SC5','B_cells_SC6')

plt_Bcells <- Meta %>% 
  mutate(ReCluster = factor(SCs)) %>%
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
  ggpubr::stat_compare_means(comparisons = my_comp, method = "t.test", label = "p.signif") +
  #ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = F, vjust = 0.5) + 
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face='bold')) +
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  
  scale_fill_manual(values=col_patient_groups) + #***
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) + #    ylab('% PBMC') + xlab('Age groups')
  ylab('% in B cells') + xlab('Clinical groups')

print(plt_Bcells)

