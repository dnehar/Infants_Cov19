
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# colors 
cols <- c('B_cells_SC0'='tomato','B_cells_SC1'='paleturquoise','B_cells_SC3'='mediumseagreen','B_cells_SC2'='cornflowerblue','B_cells_SC4'='mediumpurple')

# load meta data (PBMCs, n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')
head(Meta)
dim(Meta) 


# subset to be plotted 
subset_to_be_plotted <- c( 'B_cells_SC0','B_cells_SC1','B_cells_SC2','B_cells_SC3','B_cells_SC4')

pHC <-  paste0("pHC", seq(1:14))
pG1 <- c("pCoV1", "pCoV4", "pCoV10", "pCoV13", "pCoV14", "pCoV15", "pCoV16", "pCoV18", "pCoV19", "pCov25")
pG2 <-  c("pCoV3",  "pCoV5",  "pCoV6",  "pCoV8", "pCoV9", "pCoV12", "pCoV20", "pCov22" ,"pCov23", "pCov24", "pCov26")
pG3 <- c( "pCoV2","pCoV7", "pCoV11", "pCoV17","pCov21")
ordered_names <- c(pHC, pG1, pG2, pG3)

BP <- Meta %>% 
  
  mutate(Patient_groups = factor(Groups, levels = c("pHC", "G1", "G2", "G3"))) %>%
  mutate(ReCluster = factor(annotated_SCs)) %>%
  #mutate(ReCluster = factor(Names, levels = order_names)) %>%
  filter(ReCluster %in% subset_to_be_plotted) %>%
  
  group_by(Groups, Fig_ids, ReCluster) %>%
  #filter(Groups %in% c("HO_M",'HO_F')) %>% 
  summarise(n = n()) %>% #, Set = first(Set)
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>%
  ggplot(aes(x = Fig_ids, y = freq, fill = ReCluster, group = Groups)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=cols) + #***
  scale_x_discrete(limits=ordered_names) + #labels= labels
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle = 90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") + #    ylab('% PBMC') + xlab('Age groups') 
  ylab('% in B cells') + xlab('Individuals (n=40)')

BP 
