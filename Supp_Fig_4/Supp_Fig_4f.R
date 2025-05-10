library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)


# load meta data (PBMCs, n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# colors
col_pGroups = c('pHC'="#66bd63",'G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f")


#### ISGhi Transitional B cells (TrB; n=5228 cells)

p_tonic <- Meta %>% 
  
  filter(SCs %in% c('B_cells_SC1')) %>% #head()
  mutate(Groups = factor(Patient_groups, levels = c("pHC", "G1","G2","G3"))) %>%
  #mutate(ReCluster = factor(Final_annotations, levels = ordered_SC)) %>% #*****
  group_by(Groups, Fig_ids) %>% #SCs
  summarise(n = n()) %>% #, Set = first(Set)
  #summarise(n = n(), Age_months = first(as.numeric(Age_months)), Gender = first(as.numeric(Gender))) %>% #, Set = first(Set)
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>%
  #filter(Groups %in% c('cSLE')) %>% #head()
  arrange(desc(n)) %>% 
  mutate(Names = factor(Fig_ids, levels = unique(Fig_ids))) %>% 
  #plot the number of cells
  ggplot(aes(x=Names, y=n, fill=Groups)) + 
  geom_bar(stat="identity", color="black") + #CD4 T cells:  #697d35 #CD14 mo: #f15d64,  B cell: "#4c459c"
  ylab("Number of cells") +
  xlab("Individuals") +
  
  scale_fill_manual(values=col_pGroups) + 
  #facet_wrap(.~SCs, scales = "free_y", nrow = 2) + 
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle = 90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        panel.background = element_rect(fill = "white"),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        legend.position = "none") +
  
  ggtitle("ISGhi Tr B cells (SC1) ")  +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size=20, face="bold"))

print(p_tonic)



#### ISGhi Naive B cells (n=4985 cells)

p_tonic <- Meta %>% 
  
  filter(SCs %in% c('B_cells_SC2')) %>% #head()
  mutate(Groups = factor(Patient_groups, levels = c("pHC", "G1","G2","G3"))) %>%
  #mutate(ReCluster = factor(Final_annotations, levels = ordered_SC)) %>% #*****
  group_by(Groups, Fig_ids) %>% #SCs
  summarise(n = n()) %>% #, Set = first(Set)
  #summarise(n = n(), Age_months = first(as.numeric(Age_months)), Gender = first(as.numeric(Gender))) %>% #, Set = first(Set)
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>%
  #filter(Groups %in% c('cSLE')) %>% #head()
  arrange(desc(n)) %>% 
  mutate(Names = factor(Fig_ids, levels = unique(Fig_ids))) %>% 
  #plot the number of cells
  ggplot(aes(x=Names, y=n, fill=Groups)) + 
  geom_bar(stat="identity", color="black") + #CD4 T cells:  #697d35 #CD14 mo: #f15d64,  B cell: "#4c459c"
  ylab("Number of cells") +
  xlab("Individuals") +
  
  scale_fill_manual(values=col_pGroups) + 
  #facet_wrap(.~SCs, scales = "free_y", nrow = 2) + 
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle = 90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        panel.background = element_rect(fill = "white"),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        legend.position = "none") +
  
  ggtitle("ISGhi naive B cells (SC1) ")  +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size=20, face="bold"))

print(p_tonic)

