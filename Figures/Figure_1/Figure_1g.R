library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)
library(reshape2)

# ============================================================================
# Figure 1g: Cell Type Frequency Distribution Across Individuals by Group
# ============================================================================

# Load metadata (PBMCs, n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv') 

# Get total number of individuals
nbr_indiv <- length(unique(Meta$Fig_ids))

# ---- Order samples depending on the patient groups ----
pHC <- unique(filter(Meta, Patient_groups=='pHC') %>% select(Fig_ids))$Fig_ids
G1 <- c("pCoV1", "pCoV4", "pCoV10", "pCoV13", "pCoV14", "pCoV15", "pCoV16", "pCoV18", "pCoV19", "pCov25")
G2 <- c("pCoV3",  "pCoV5",  "pCoV6",  "pCoV8", "pCoV9", "pCoV12", "pCoV20", "pCov22", "pCov23", "pCov24", "pCov26")
G3 <- c("pCoV2","pCoV7", "pCoV11", "pCoV17","pCov21")
ordered_names <- c(pHC, G1, G2, G3)

# ---- Define color schemes ----
# Group-level colors for x-axis labels
group_colors <- c(rep("#66bd63", length(pHC)),      # pHC (green)
                  rep("#d8daeb", length(G1)),       # G1 (light purple)
                  rep("#9e9ac8", length(G2)),       # G2 (medium purple)
                  rep("#54278f", length(G3)))       # G3 (dark purple)

# Cell type colors
col_simple_clustering <- c(
  'CD8_Tcells' = '#fc8d59',
  'CD4_Tcells' = '#9e9ac8',
  'B_cells' = '#96daf7',
  'NK' = '#fed976',
  'HSCs' = '#b0479a',
  'Eryth' = '#35978f',
  'pDC' = '#193a1c',
  'Mgk' = '#8c510a',
  'cDC' = '#e31a1c',
  'CD14_mono' = '#f6a2a7',
  'CD16_mono' = '#f9d3d7',
  'PCs' = '#08306b'
)

# Order clusters for stacking
ordered_clusters <- c('CD14_mono', 'CD16_mono', 'cDC', 'pDC', 'Mgk', 'NK', 'CD8_Tcells',
                      'CD4_Tcells', 'B_cells', 'PCs', 'HSCs', 'Eryth')

# ---- Calculate cell type proportions per individual ----
y <- prop.table(x = table(as.character(Meta$simple_clustering), Meta$Fig_ids), margin = 2)
df1 <- reshape2::melt(y)
colnames(df1) <- c("clusters", "Indiv", "Percentage")
df1 <- df1 %>% mutate(clusters = factor(clusters, levels = ordered_clusters))

# Display data summary
cat("Data dimensions:", dim(df1), "\n")
print(head(df1))

# ---- Create stacked bar plot ----
p_Ind <- ggplot(data = df1, aes(x = as.character(Indiv), y = Percentage, fill = clusters)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  scale_fill_manual(values = col_simple_clustering) +
  scale_x_discrete(limits = ordered_names) +
  theme_minimal() +
  theme(
    # X-axis formatting
    axis.text.x = element_text(colour = group_colors, face = "bold", size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(fill = NA, color = 'black', size = 1),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold")
  ) +
  xlab("Individuals (n=40)") +
  ylab("Cell frequency")

print(p_Ind)

# ---- Export figure ----
ggsave("./PANELS/Barplot_pCov_cluster_perIndiv_11272024.pdf", 
       plot = p_Ind,
       width = 3.4, height = 1.4, units = "in", scale = 3,
       dpi = 300)

cat("Figure saved successfully!\n")
