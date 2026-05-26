# Load required libraries for data manipulation and visualization
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Define color scheme for clinical groups (purple gradient for patient groups, green for healthy controls)
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Load and inspect metadata for PBMC cells (n=203,402 cells)
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')
head(Meta)
dim(Meta) 

# Define the clinical groups to compare and generate all pairwise comparisons
clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups, 2, FUN = list, simplify = T)

# Specify which dendritic cell subtypes to include in the plot
subset_to_be_plotted <- c( 'cDC_SC0','cDC_SC1','cDC_SC2','cDC_SC3','cDC_SC4','cDC_SC5','cDC_SC6')

# Create boxplot with jittered points showing DC subset frequencies across clinical groups
plt_clinical <- Meta %>% 
  # Convert cluster annotations to factor
  mutate(ReCluster = factor(annotated_SCs)) %>%
  # Convert clinical groups to factor with specified order
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%
  # Filter to only include specified DC subtypes
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Group by clinical group, figure ID, and cluster, then count cells
  group_by(Groups, Fig_ids, ReCluster) %>%
  summarise(n = n()) %>% #, Set = first(Set)
  # Calculate frequency as percentage of total cells per group
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>% 
  # Create ggplot visualization with boxplots and overlaid points
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  # Draw boxplots (outliers hidden for clarity)
  geom_boxplot(outlier.shape = NA) +
  # Add jittered individual data points on top of boxplot
  geom_jitter(size = 0.2) +
  theme_bw()  +  #THEME +
  # Perform pairwise t-tests and display significance markers (*, **, ***)
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  # Remove legend and format facet labels
  theme(legend.position = "none", 
        strip.text = element_text(size = 14, face='bold')) +
  # Create separate facets for each DC subtype in one row with independent y-axes
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) + 
  # Apply color scheme to fill boxplots
  scale_fill_manual(values=col_pGroups) + #**
  # Format axis text and labels
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =0),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) + #    ylab('% PBMC') + xlab('Age groups')
  # Set axis labels
  ylab('% in DCs') + xlab('Clinical groups')

# Display the final plot
print(plt_clinical)
