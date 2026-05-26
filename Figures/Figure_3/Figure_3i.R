# Load required libraries for data manipulation and visualization
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Define color palette for patient groups: G1 (light purple), G2 (medium purple), G3 (dark purple), pHC (green)
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Load metadata containing CD4+ T regulatory cell (Treg) subclusters
Meta <- read.csv('./meta_subclusters/pCoV_CD4_Tregs_bbknn_02042025.csv')
head(Meta)
dim(Meta)

# Create the clinical groups boxplot by:
plt_clinical <- Meta %>% 
  # Convert subclusters to factor for categorical grouping
  mutate(ReCluster = factor(SCs)) %>%
  # Convert patient groups to ordered factor using predefined levels
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Group data by clinical groups, figure IDs, and reclusters
  group_by(Groups, Fig_ids, ReCluster) %>%
  # Count cells and calculate frequency (percentage) within each group
  summarise(n = n()) %>%
  mutate(freq = n / sum(n) *100) %>%
  # Ungroup data for plotting
  ungroup() %>%
  as.data.frame() %>% 
  # Initialize ggplot with groups on x-axis, frequency on y-axis, colored by groups
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  # Add boxplots without outlier markers
  geom_boxplot(outlier.shape = NA) +
  # Overlay individual data points with transparency
  geom_jitter(size = 0.2) +
  theme_bw() +
  # Perform statistical comparison between groups using t-tests, showing significance symbols
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend and format facet titles
  theme(legend.position = "none", strip.text = element_text(size = 14, face='bold')) +
  # Create separate panels (facets) for each reCluster with independent y-axis scales
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  # Apply custom colors to the groups
  scale_fill_manual(values=col_pGroups) +
  # Format axis text and labels for readability
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle=0),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) +
  # Set axis labels
  ylab('% in Tregs') + xlab('Clinical groups')

# Display the plot
print(plt_clinical)
