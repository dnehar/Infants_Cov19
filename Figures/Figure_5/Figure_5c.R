# Load required libraries
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)


# Load metadata containing PBMC cell information (n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# Define color palette for clinical groups
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Define clinical group order and generate all pairwise comparisons for statistical testing
clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups, 2, FUN = list, simplify = T)

# Define B cell subclusters to include in the plot
subset_to_be_plotted <- c('B_cells_SC0','B_cells_SC1','B_cells_SC2','B_cells_SC3','B_cells_SC4')

# Create boxplot with jittered points showing B cell subcluster frequencies across clinical groups
plt_clinical <- Meta %>% 
  # Convert SCs to factor with specified levels to ensure correct order
  mutate(ReCluster = factor(SCs, levels = subset_to_be_plotted)) %>%
  # Convert clinical groups to factor with specified levels
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Filter to only include the specified B cell subclusters
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Group by clinical group, sample ID, and subcluster to count cells
  group_by(Groups, Fig_ids, ReCluster) %>%
  summarise(n = n()) %>%
  # Calculate the percentage frequency of each subcluster within B cells for each sample
  mutate(freq = n / sum(n) * 100) %>%
  ungroup() %>%
  as.data.frame() %>% 
  # Create the ggplot with boxplot and overlay jittered points
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.2) +
  theme_bw() +
  # Add statistical significance testing (t-tests) between all pairwise combinations
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend and format facet labels
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face = 'bold')) +
  # Create separate panels for each B cell subcluster with independent y-axis scales
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  # Apply custom color scheme to each clinical group
  scale_fill_manual(values = col_pGroups) +
  # Format axis labels and titles for clarity
  theme(axis.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 16, angle = 0),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18)) +
  # Set axis and plot labels
  ylab('% in B cells') + 
  xlab('Clinical groups')

# Display the final plot
print(plt_clinical)
