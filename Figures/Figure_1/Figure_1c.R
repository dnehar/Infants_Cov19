# Load required libraries for data manipulation and visualization
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Load metadata file containing PBMC cell annotations (203,402 cells total)
Meta <- read.csv('Meta_pCoV40_03112025_small.csv') 

# Define a color palette for cell type visualization
# Maps each cell type/cluster name to a specific hex color code
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

# Create UMAP visualization colored by cell type clusters
p_umap1 <- ggplot(Meta) +
  # Plot each cell as a point with UMAP coordinates, colored by cluster assignment
  geom_point(aes(x=X_umap1, y=X_umap2, color=factor(simple_clustering)), 
             size=0.2) +
  # Apply the custom color palette to the plot
  scale_color_manual(values=col_simple_clustering, 'Sublusters')+
  # Label the axes
  xlab('UMAP_1') +
  ylab('UMAP_2') + 
  # Remove default ggplot theme elements for a clean look
  theme_void() +
  # Format the plot title (bold, larger font)
  theme(plot.title = element_text(size = 20, face = "bold", vjust = 0.03)) +
  # Format legend text (smaller, bold)
  theme(legend.text=element_text(size=7, face = "bold"))+
  # Increase legend point size for visibility
  guides(color = guide_legend(override.aes = list(size=7))) +
  # Add plot title
  labs(title = "annotated clusters") 

# Display the plot
p_umap1
