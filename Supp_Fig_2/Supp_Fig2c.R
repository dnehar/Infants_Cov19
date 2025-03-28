library(Seurat)
library(dplyr)
library(cowplot)
library(tidyr)
library(Nebulosa)
library(anndata)
   
# load anndta object object
adata <- read_h5ad("./H5AD_obj/CD14_mono_cleaned_pCov40_11232021.h5ad")
adata

# Generate seurat object #########
    # raw counts 
   raw_exprs <- t(as.matrix(adata$raw$X))
   colnames(raw_exprs) <- adata$raw$obs_names #***
   rownames(raw_exprs) <- adata$raw$var_names
   
   
   # Create the Seurat object
   seurat <- CreateSeuratObject(raw_exprs)
   #seurat <- CreateSeuratObject(exprs)
   
   # Set the expression assay
   #seurat <- SetAssayData(seurat, "data", exprs)
   
   seurat <- SetAssayData(seurat, "data", raw_exprs)
   
   # Add observation metadata
   seurat <- AddMetaData(seurat, adata$obs)
   
   # Add fetaure metadata
   seurat[["RNA"]][["n_cells"]] <- adata$var["n_cells"]
   
   # Add UMAP embedding
   embedding <- data.frame(adata$obsm["X_umap"])
   rownames(embedding) <- adata$obs_names #$to_list()
   colnames(embedding) <- c("umap_1", "umap_2")#) #,"umap_3"
   head(embedding)
   seurat[["umap"]] <- CreateDimReducObject(as.matrix(embedding), key = "umap_")
   
  # Define identity 
   seurat <- SetIdent(seurat, value = "SCs")
   DimPlot(seurat)
 
 # plot genes using Nebulosa ##########
 
 p_CD16_mo <- plot_density(seurat, 
                             features =  c('FCGR3A','S100A12','IFI27','MX1','IFI44L','IFI44'),
                             size = 0.5, pal='inferno',
                             slot = 'counts') &  theme_void()  & theme(plot.title = element_text(size = 20, face = 'bold',hjust = 0.5,family = 'Helvetica')
                             )  
   print(p_CD16_mo)
   
   # save figure 
   ggsave("./PANELS/CD16_mo_nebulosa_pCoV_12022024.pdf", p_CD16_mo,
          width=6.8, height=3.8,  units="in", scale=1)
   
