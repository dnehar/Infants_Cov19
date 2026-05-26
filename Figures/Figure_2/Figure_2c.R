# Load required libraries for data manipulation, plotting, and statistical comparison
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)


# Load metadata containing PBMC cell information (n=203,402 cells)
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# Define color palette for patient groups
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Define clinical group order and generate all pairwise comparisons
clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups, 2, FUN = list, simplify = T)

# Define monocyte subclusters to visualize
subset_to_be_plotted <- c( 'CD14_mono_SC0','CD14_mono_SC1','CD14_mono_SC2', 'CD14_mono_SC3')

# Create boxplot comparing frequencies of monocyte subclusters across clinical groups
plt_clinical <- Meta %>% 
  # Convert cluster and group annotations to factors
  mutate(ReCluster = factor(annotated_SCs)) %>%
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Filter to only include specified monocyte subclusters
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Count cells per group and cluster, grouping by patient ID (Fig_ids)
  group_by(Groups, Fig_ids, ReCluster) %>%
  summarise(n = n()) %>%
  # Calculate frequency as percentage within each cluster
  mutate(freq = n / sum(n) * 100) %>%
  ungroup() %>%
  as.data.frame() %>% 
  # Create ggplot visualization
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  # Add boxplot layer with no outliers shown separately
  geom_boxplot(outlier.shape = NA) +
  # Overlay individual data points
  geom_jitter(size = 0.2) +
  theme_bw() +
  # Add statistical comparison (t-test) with significance indicators between groups
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend and format facet labels
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face='bold')) +
  # Create separate panels for each monocyte subcluster
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  # Apply custom color palette
  scale_fill_manual(values=col_pGroups) +
  # Format axis labels and titles
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =0),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) +
  # Add axis labels
  ylab('% in CD14 mo') + xlab('Clinical groups')

# Display the plot
print(plt_clinical)
