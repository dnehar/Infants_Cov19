
  library(dplyr)
  library(cowplot)
  library(ggpubr)
  library(tidyr)

# load meta data (PBMCs, n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv') 


# colors 
col_simple_clustering = c('CD8_Tcells' = '#fc8d59',
                          'CD4_Tcells' = '#9e9ac8',
                          'B_cells'='#96daf7',
                          'NK' = '#fed976',
                          'HSCs' = '#b0479a',
                          'Eryth' = '#35978f',
                          'pDC' = '#193a1c',
                          #'Plasma' = '#232323',
                          'Mgk' = '#8c510a',
                          'cDC' = '#e31a1c',
                          'CD14_mono' = '#f6a2a7',
                          'CD16_mono' = '#f9d3d7',
                          'PCs' = '#08306b')



p_umap1 <- ggplot(Meta) +
  geom_point(aes(x=X_umap1, y=X_umap2,  color=factor(simple_clustering)), 
             size=0.2) +
  scale_color_manual(values=col_simple_clustering,'Sublusters')+
  xlab('UMAP_1') +
  ylab('UMAP_2') + 
  theme_void() +
  theme(plot.title = element_text(size = 20, face = "bold", vjust = 0.03)) +
  theme(legend.text=element_text(size=7, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size=7))) +
  labs(title = "annotated clusters") 

p_umap1
