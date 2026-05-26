# Load required libraries for data manipulation, visualization, and statistical comparison
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Load metadata containing PBMC cell information (203,402 cells total)
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# Define color palette for the four clinical groups (pHC=healthy control, G1-G3=patient groups)
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Define the order of clinical groups and generate all pairwise comparisons for statistical testing
clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups, 2, FUN = list, simplify = T)

# Select the three plasma cell subclusters to plot
subset_to_be_plotted <- c('PCs_SC0', 'PCs_SC1', 'PCs_SC2')

# Create the plot: start with metadata and apply a series of transformations
plt_clinical <- Meta %>% 
  # Convert subcluster column to an ordered factor
  mutate(ReCluster = factor(SCs, levels = subset_to_be_plotted)) %>%
  # Convert patient groups column to an ordered factor
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Keep only rows matching the selected subclusters
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Group data by clinical group, figure ID, and subcluster to count cells
  group_by(Groups, Fig_ids, ReCluster) %>%
  # Count cells in each group
  summarise(n = n()) %>%
  # Calculate frequency as percentage within each group
  mutate(freq = n / sum(n) * 100) %>%
  # Ungroup the data
  ungroup() %>%
  as.data.frame() %>% 
  
  # Initialize ggplot with clinical groups on x-axis and frequency on y-axis
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  # Add boxplot layer (hide outliers as they're shown by jitter)
  geom_boxplot(outlier.shape = NA) +
  # Overlay individual data points with slight transparency
  geom_jitter(size = 0.2) +
  # Use black-and-white theme
  theme_bw() +
  # Add statistical comparison (t-test) between all group pairs with significance stars
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend and bold facet titles
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face = 'bold')) +
  # Create separate facets for each subcluster (1 row, free y-axis scaling)
  facet_wrap(. ~ ReCluster, scales = "free_y", nrow = 1) +
  # Apply custom colors to each clinical group
  scale_fill_manual(values = col_pGroups) +
  # Customize axis labels and text sizes
  theme(axis.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 16, angle = 0),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18)) +
  # Add axis labels
  ylab('% in plasma cells (PCs)') + xlab('Clinical groups')

# Display the plot
print(plt_clinical)
