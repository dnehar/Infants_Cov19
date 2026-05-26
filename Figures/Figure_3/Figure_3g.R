# Load required libraries
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Define color palette for patient groups
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Load metadata file containing clinical group and cluster information
Meta <- read.csv('./meta_subclusters/pCoV_CD4_Tmem_02042025.csv')

# Create boxplot with statistical comparison
plt_clinical <- Meta %>% 
  # Convert subcluster column to factor
  mutate(ReCluster = factor(SCs)) %>%
  # Convert patient groups to factor with specified order
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Group by clinical groups, figure IDs, and reclusters
  group_by(Groups, Fig_ids, ReCluster) %>%
  # Count cells in each group
  summarise(n = n()) %>%
  # Calculate frequency as percentage of total cells
  mutate(freq = n / sum(n) *100) %>%
  # Ungroup and convert to data frame
  ungroup() %>%
  as.data.frame() %>% 
  # Initialize ggplot with groups on x-axis, frequency on y-axis
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  # Add boxplot layer with hidden outliers
  geom_boxplot(outlier.shape = NA) +
  # Overlay individual data points with jitter to avoid overplotting
  geom_jitter(size = 0.2) +
  # Apply black and white theme
  theme_bw()  +
  # Add statistical comparison between groups using t-test with significance symbols
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend and format facet labels
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face='bold')) +
  # Create separate panels for each cluster with free y-axis scaling
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  # Apply custom color palette to fill
  scale_fill_manual(values=col_pGroups) +
  # Format axis labels and titles (size and bold)
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =0),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) +
  # Set axis labels
  ylab('% in Tmem') + xlab('Clinical groups')

# Display the plot
print(plt_clinical)
