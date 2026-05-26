# Load required libraries for data manipulation and visualization
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)


# Load metadata containing cell annotations and clinical information (203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# Define color palette for clinical groups (G1, G2, G3 = patient groups; pHC = healthy controls)
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Specify the three CD16 monocyte subclusters to plot
subset_to_be_plotted <- c( 'CD16_mono_SC0','CD16_mono_SC1','CD16_mono_SC2')

# Create boxplot visualization comparing CD16 monocyte subcluster frequencies across clinical groups
plt_clinical <- Meta %>% 
  # Convert annotated subclusters to factor for plotting
  mutate(ReCluster = factor(annotated_SCs)) %>%
  # Convert patient groups to ordered factor based on predefined clinical groups
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Filter to include only the three CD16 monocyte subclusters
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Count cells per clinical group, patient, and subcluster
  group_by(Groups, Fig_ids, ReCluster) %>%
  summarise(n = n()) %>%
  # Calculate frequency (%) of each subcluster within each group
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>% 
  # Create ggplot with clinical groups on x-axis, frequency on y-axis
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  # Add boxplot layer with outliers hidden
  geom_boxplot(outlier.shape = NA) +
  # Overlay individual data points with slight jitter to show variability
  geom_jitter(size = 0.2) +
  # Apply clean black and white theme
  theme_bw()  +
  # Add statistical comparison between groups using t-tests with significance symbols
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend since clinical groups are already labeled on x-axis
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face='bold')) +
  # Create separate subplots for each CD16 monocyte subcluster with independent y-axes
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  # Apply custom color palette to fill boxplots
  scale_fill_manual(values=col_pGroups) +
  # Format axis labels and titles with larger font sizes and bold titles
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =0),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) +
  # Set axis labels
  ylab('% in CD16 mo') + xlab('Clinical groups')

# Display the plot
print(plt_clinical)
