library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Load metadata containing PBMC cell data (n=203,402 cells)
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# Define color palette for clinical groups
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Set the order of clinical groups for plotting
clinical_groups <- c("pHC", "G1", "G2", "G3")

# Generate all pairwise comparisons between clinical groups for statistical testing
my_comparisons <- combn(clinical_groups, 2, FUN = list, simplify = T)

# Define the CD8 T cell subclusters to be visualized
subset_to_be_plotted <- c('CD8_NAIVE', 'CD8_ISGhi', 'CD8_GzK', 'CD8_TEMRA', 'CD8_Prolif')

# Create the plot
plt_clinical <- Meta %>% 
  # Convert annotated_SCs to factor with specified order
  mutate(ReCluster = factor(annotated_SCs, levels = subset_to_be_plotted)) %>%
  # Convert Patient_groups to factor with specified clinical group order
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Filter to keep only the CD8 subclusters of interest
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Group by clinical group, patient ID, and CD8 subcluster
  group_by(Groups, Fig_ids, ReCluster) %>%
  # Count cells in each group
  summarise(n = n()) %>%
  # Calculate frequency as percentage within each clinical group
  mutate(freq = n / sum(n) * 100) %>%
  # Ungroup for ggplot
  ungroup() %>%
  as.data.frame() %>% 
  # Create boxplot with individual data points
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  # Add boxplot layer without outliers (shown as individual points instead)
  geom_boxplot(outlier.shape = NA) +
  # Overlay individual data points with jitter
  geom_jitter(size = 0.2) +
  theme_bw() +
  # Add statistical significance testing (t-test) between all pairwise comparisons
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend and format facet labels
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face = 'bold')) +
  # Create separate subplots for each CD8 subcluster in a single row
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  # Apply custom color palette to fill aesthetic
  scale_fill_manual(values = col_pGroups) +
  # Customize axis text and labels
  theme(axis.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 16, angle = 0),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18)) +
  # Label axes
  ylab('% in CD8 T cells') +
  xlab('Clinical groups')

# Display the plot
print(plt_clinical)
