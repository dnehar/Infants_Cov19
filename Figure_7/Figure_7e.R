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



 # colors  
  col_SCs <- c('B_cells_SC0'='tomato',
              'B_cells_SC1'='mediumseagreen',
              'B_cells_SC2'='cornflowerblue',
              'B_cells_SC3'='paleturquoise',
              'B_cells_SC4'='mediumpurple',
              'B_cells_SC5'='goldenrod',
              'B_cells_SC6'='lightgreen')
              

# B cell subsets 
subset_to_be_plotted <- c("B_cells_SC0", paste0("B_cells_SC",seq(1:5)))
  
# order names 
 pHC <-  paste0("pHC", seq(1:14))
  pG1 <- c("pCoV1", "pCoV4", "pCoV10", "pCoV13", "pCoV14", "pCoV15", "pCoV16", "pCoV18", "pCoV19", "pCov25")
  pG2 <-  c("pCoV3",  "pCoV5",  "pCoV6",  "pCoV8", "pCoV9", "pCoV12", "pCoV20", "pCov22" ,"pCov23", "pCov24", "pCov26")
  pG3 <- c( "pCoV2","pCoV7", "pCoV11", "pCoV17","pCov21")
  aHC <-  c(paste0("aHC", seq(1:3)), paste0("aHC", seq(7, 14, by=1))) 
  aG1 <- c("aCoV8","aCoV9","aCoV20")
  aG2 <- c("aCoV3", "aCoV21", "aCoV25", "aCoV28"  )
  aG3 <- c("aCoV1", "aCoV2" ,"aCoV4",  "aCoV5",  "aCoV6",   "aCoV7", "aCoV10", "aCoV11", "aCoV12", "aCoV13",
           "aCoV14", "aCoV15", "aCoV16", "aCoV17", "aCoV18", "aCoV19",  "aCoV22", "aCoV23", "aCoV24", "aCoV26", "aCoV27") 
                                  
  ordered_names <- c(pHC, pG1, pG2, pG3, aHC, aG1, aG2, aG3)
  length(ordered_names)
                                   
  
  BP <- Meta %>% 
    
    mutate(Patient_groups = factor(Groups, levels = clinical_groups)) %>%
    mutate(ReCluster = factor(SCs)) %>%
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
    scale_fill_manual(values=col_SCs) + #***
    scale_x_discrete(limits=ordered_names) + #labels= labels
    theme(axis.text.y=element_text(size=16), 
          axis.text.x=element_text(size=16, angle = 90),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          legend.position = "none") + #    ylab('% PBMC') + xlab('Age groups')
    
    ylab('% in CD8 T cells') + xlab('Individuals (n=40)')
  
  print(BP)
