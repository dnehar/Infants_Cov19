# Load required libraries for data manipulation and visualization
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Define color palette for clinical groups
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Load metadata containing PBMC cell information (203,402 cells total)
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')
head(Meta)
dim(Meta) 

# Define clinical group order and generate all pairwise comparisons for statistical testing
clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups,2, FUN = list, simplify = T)

# Define CD4 T cell subclusters to include in the plot
subset_to_be_plotted <- c( 'CD4_SC0','CD4_SC1','CD4_SC2','CD4_SC3')

# Create the main plot: boxplot with jitter points showing cell frequency by group
plt_clinical <- Meta %>% 
  # Convert cell cluster annotations to factor
  mutate(ReCluster = factor(annotated_SCs)) %>%
  # Convert patient groups to ordered factor
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Keep only CD4 subclusters of interest
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Count cells per group and cluster combination
  group_by(Groups, Fig_ids, ReCluster) %>%
  summarise(n = n()) %>%
  # Calculate frequency as percentage
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>% 
  # Initialize ggplot with groups on x-axis and frequency on y-axis
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  # Add boxplot layer (hide outliers for clarity)
  geom_boxplot(outlier.shape = NA) +
  # Overlay individual data points with slight jitter
  geom_jitter(size = 0.2) +
  # Use clean black and white theme
  theme_bw()  +
  # Add statistical significance testing using t-tests with pairwise comparisons
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend and bold facet labels
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face='bold')) +
  # Create separate panels for each CD4 subcluster
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  # Apply custom color palette
  scale_fill_manual(values=col_pGroups) +
  # Customize text sizes and styling for axes
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =0),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) +
  # Add axis labels
  ylab('% in CD4 T cells') + xlab('Clinical groups')

# Display the plot
print(plt_clinical)
