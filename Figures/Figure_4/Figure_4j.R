# Load required libraries
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Load metadata for PBMC cells (n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# Define color palette for clinical groups
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Define clinical groups and generate all pairwise comparisons for statistical testing
clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups, 2, FUN = list, simplify = T)

# Define NK cell subclusters to be analyzed and plotted
subset_to_be_plotted <- c('NK_SC0', 'NK_SC1', 'NK_SC2', 'NK_SC3')

# Create main plot showing frequency of NK cell subclusters across clinical groups
plt_clinical <- Meta %>% 
  # Convert subcluster and group columns to factors with specified order
  mutate(ReCluster = factor(SCs, levels = subset_to_be_plotted)) %>%
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Filter to include only the NK cell subclusters of interest
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Count cells in each subcluster per patient per group
  group_by(Groups, Fig_ids, ReCluster) %>%
  summarise(n = n()) %>%
  # Calculate frequency as percentage within each NK cell subcluster
  mutate(freq = n / sum(n) * 100) %>%
  ungroup() %>%
  as.data.frame() %>%
  # Create ggplot boxplot with individual points (jitter)
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.2) +
  theme_bw() +
  # Add statistical comparisons (t-tests with significance labels)
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend and customize strip text
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face = 'bold')) +
  # Create separate facets for each NK cell subcluster
  facet_wrap(. ~ ReCluster, scales = "free_y", nrow = 1) +
  # Apply custom color palette
  scale_fill_manual(values = col_pGroups) +
  # Customize axis labels and text sizes
  theme(axis.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 16, angle = 0),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18)) +
  # Add axis labels
  ylab('% in NK cells') + 
  xlab('Clinical groups')

# Display the plot
print(plt_clinical)
