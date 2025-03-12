library(matrixTests)
library(tidyverse)
library(ComplexHeatmap)

npx <- read.csv("olink_raw_data.csv")
group <- sapply(npx$Sample,str_split_i, pattern = "_",i=1) %>% unname
group <- factor(group, c("Healthy","G1","G2","G3"))
row.names(npx) <- npx$Sample
npx <- t(npx[,-1])

anova.res <- row_oneway_equalvar(npx,group)
anova.res$padj <- p.adjust(anova.res$pvalue)
sig.cytokine <- anova.res %>% filter(padj<0.1) %>% row.names

mat <- npx[sig.cytokine,]
mat <- mat - rowMedians(mat[,1:20])

ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill =c( "#78C679", "#D8DAEB", "#9E9AC8", "#54278F") ),
                                         labels = c("pHC", "G1", "G2","G3"), 
                                         labels_gp = gpar(col = "white", fontsize = 10)))

Heatmap(as.matrix(mat), 
        column_title = "Cytokine", 
        col=circlize::colorRamp2(c(-3,0,3), c("#FF00FF","#000000", "#FFFF00")),
        na_col = "white",
        top_annotation = ha,
        row_dend_side = "left",
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        cluster_columns = F,
        show_row_names = T,
        show_column_names = F,
        column_split = group,
        name = "NPX"
)
