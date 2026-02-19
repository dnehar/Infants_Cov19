library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)


# load meta data (PBMCs, n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv') 

nbr_indiv <- length(unique(SC_Meta$Fig_ids))
  
  #---- Order samples depending on the patient groups 
pHC <- unique(filter(SC_Meta, Patient_groups=='pHC') %>% select(Fig_ids))$Fig_ids
G1 <- c("pCoV1", "pCoV4", "pCoV10", "pCoV13", "pCoV14", "pCoV15", "pCoV16", "pCoV18", "pCoV19", "pCov25")
G2 <-  c("pCoV3",  "pCoV5",  "pCoV6",  "pCoV8", "pCoV9", "pCoV12", "pCoV20", "pCov22" ,"pCov23", "pCov24", "pCov26")
G3 <- c( "pCoV2","pCoV7", "pCoV11", "pCoV17","pCov21")
ordered_names <- c(pHC, G1, G2, G3)
  
#- colors
  a <- c(rep("#66bd63",length(pHC)),#14
         rep("#d8daeb",length(G1)), #10
         rep("#9e9ac8",length(G2)),#11
         rep("#54278f",length(G3)))# 5
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
  
  # order clusters
  ordered_clusters <- c('CD14_mono','CD16_mono','cDC','pDC','Mgk','NK','CD8_Tcells',
                        'CD4_Tcells','B_cells','PCs','HSCs','Eryth')
  
  
  #head(SC_Meta)
  #- SC_ids column for for Barplot per donor 
  y=prop.table(x=table(as.character(Meta$simple_clustering), Meta$Fig_ids), margin=2)
  df1 <- melt(y)
  colnames(df1) <- c("clusters", "Indiv", "Percentage")
  df1 %>% mutate(clusters = factor(clusters, levels = ordered_clusters)) -> df1
  
  head(df1)
  dim(df1)
  
  p_Ind <- ggplot(data=df1, aes(x=as.character(Indiv), y=Percentage, fill=clusters)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values=col_simple_clustering) +
    #scale_x_discrete(limits=c(as.character(seq(1:(nbr_indiv))))) + #labels= labels
    scale_x_discrete(limits=ordered_names) + #labels= labels
    #THEME +
    theme(axis.text.x = element_text(colour = a)) + #### x axis color per group
    ylab("Cell frequency") +
    #guides(col = guide_legend(ncol = 2)) +
    theme(axis.text.y=element_text(face="bold",size=14), 
          axis.text.x=element_text(face="bold",size=14, angle = 90),
          axis.title.x = element_text(face="bold", size=16),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(hjust = 0.5,face='bold',size=14),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, color = 'black', size=1),
          panel.background = element_rect(fill = "white", colour = NA),
          legend.text = element_text(size=9)) +
    #theme(legend.position="none") +
    
    xlab("Individuals (n=40)") #+
  
  #ggtitle("Individuals") +
  #ggtitle(paste(Subcluster,"SC abundance across individuals"))
  print(p_Ind)
  
  
  getwd()
  ggsave("./PANELS/Barplot_pCov_cluster_perIndiv_11272024.pdf", p_Ind,
         width=3.4, height=1.4,  units="in", scale=3)
