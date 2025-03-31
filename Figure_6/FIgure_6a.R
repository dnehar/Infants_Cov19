library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)
library(pheatmap)


#- Sample information 
PhenoData <- read.table('pCov40_SampleInfo_03112025.txt', header = T, row.names = 2)
PhenoData %>% select(Groups, Patient_groups, Gender) -> Pheno  # Steroids, Vasoactive_drugs,

#- Meta data (PBMCs) 
meta <- read.csv("./Meta/Meta_pCoV40_cleaned_12062021.csv", row.names = 1)
head(meta)

y=prop.table(x=table(meta$Fig_ids, as.character(meta$annotated_SCs)), margin=2)
head(y)

# colors 
col <- colorRampPalette(c("#FF00FF","#000000", "#FFFF00"))(n=100)
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")
ann_colors=list(Groups=c(pCov="#9e9ac8", pHC="#66bd63"),
                #Batches= col_Batches,
                Patient_groups= col_pGroups)


# plot heatmap
p_HM <- pheatmap(t(y), 
                 scale = 'row', 
                 color = col, 
                 #breaks=bk,
                 cellheight = 12,  
                 cellwidth = 8, 
                 border_color = F, 
                 annotation_colors =  ann_colors, 
                 annotation_col = Pheno, 
                 #annotation_row = ann_row,
                 clustering_method = "complete", #ward.D2 
                 treeheight_row = 5, 
                 treeheight_col = 5, 
                 fontsize_row = 8,
                 fontsize_col = 8,
                 cutree_cols = 2,
                 cutree_rows = 2)
print (p_HM)
# save 
ggsave(paste0("./phearmap_SCs.pdf"), 
       p_HM,  
       width=3.5, height=3.2,  units="in", scale=3)

