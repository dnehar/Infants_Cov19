# Load required libraries for data manipulation and visualization
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Load metadata file containing PBMC cell data (203,402 cells)
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# Define color palette for the clinical groups
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Set the order of clinical groups and generate all pairwise comparisons
clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups, 2, FUN = list, simplify = T)

# Define B cell subclusters to include in the plot
subset_to_be_plotted <- c( 'NK_SC0','NK_SC1','NK_SC2','NK_SC3')

# Create plot: start with metadata and perform data preparation
plt_clinical <- Meta %>% 
  # Convert subclusters to ordered factors
  mutate(ReCluster = factor(SCs, levels = subset_to_be_plotted)) %>%
  # Convert patient groups to ordered factors
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Keep only the B cell subclusters of interest
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Group by clinical group, sample ID, and subcluster
  group_by(Groups, Fig_ids, ReCluster) %>%
  # Count cells in each group
  summarise(n = n()) %>%
  # Calculate frequency (percentage within each B cell subcluster)
  mutate(freq = n / sum(n) * 100) %>%
  # Remove grouping structure
  ungroup() %>%
  as.data.frame() %>% 
  
  # Initialize ggplot with clinical groups on x-axis and frequency on y-axis
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  # Add boxplot layer (hide outliers for clarity)
  geom_boxplot(outlier.shape = NA) +
  # Add individual data points with slight transparency
  geom_jitter(size = 0.2) +
  # Use black and white theme
  theme_bw() +
  # Add statistical comparisons (t-tests between all group pairs) with significance symbols
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend and format facet labels
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face = 'bold')) +
  # Create separate plots for each B cell subcluster (1 row, multiple columns)
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) +
  # Apply custom colors to the groups
  scale_fill_manual(values = col_pGroups) +
  # Format axis labels and titles
  theme(axis.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 16, angle = 0),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18)) +
  # Set axis labels
  ylab('% in NK cells') + xlab('Clinical groups')

# Display the plot
print(plt_clinical)
