
# Load required libraries for data manipulation and visualization
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Define color palette for the four clinical groups
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Load metadata for PBMC dataset (203,402 cells)
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')
head(Meta)
dim(Meta) 

# Define clinical groups and generate all pairwise comparisons for statistical testing
clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups,2, FUN = list, simplify = T)

# Define CD8 T cell subsets to include in the plot
subset_to_be_plotted <- c( 'CD8_NAIVE','CD8_ISGhi','CD8_GzK','CD8_TEMRA','Prolif_Cytotox', 'Prolif_GzK','Prolif_Naive')

# Create the main plot with statistical comparisons
plt_clinical <- Meta %>% 
  # Convert cell clusters to ordered factors
  mutate(ReCluster = factor(SCs, levels =subset_to_be_plotted)) %>%
  # Convert patient groups to ordered factors
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Keep only rows with the specified cell types
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Count cells per group/patient/cluster
  group_by(Groups, Fig_ids, ReCluster) %>%
  summarise(n = n()) %>%
  # Calculate frequency as percentage within each group
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>% 
  # Build ggplot with clinical groups on x-axis, frequency on y-axis
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  # Add boxplots without outlier points
  geom_boxplot(outlier.shape = NA) +
  # Overlay individual data points
  geom_jitter(size = 0.2) +
  theme_bw() +
  # Add statistical significance labels (t-test p-values) between group comparisons
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Hide legend and format facet labels
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face='bold')) +
  # Create separate subplots for each cell type with independent y-axes
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) +
  # Apply custom color palette
  scale_fill_manual(values=col_pGroups) +
  # Format axis text and labels
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =0),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) +
  # Add axis labels
  ylab('% in CD8 T cells') + xlab('Clinical groups')

# Display the final plot
plt_clinical
