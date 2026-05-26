library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Define colors for each clinical group
col_pGroups = c('G1'="#d8daeb", 'G2'="#9e9ac8", 'G3'="#54278f",'pHC'="#66bd63")

# Load cell-level metadata (PBMCs, n=203,402 cells)
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')
head(Meta)
dim(Meta) 

# Define the ordered clinical groups and all pairwise comparisons for statistical testing
clinical_groups <- c("pHC", "G1", "G2", "G3")
my_comparisons <- combn(clinical_groups,2, FUN = list, simplify = T)

# Specify the CD4 T cell subclusters to include in the plot
subset_to_be_plotted <- c( 'CD4_SC0','CD4_SC1','CD4_SC2','CD4_SC3')

# Build the boxplot: compute per-patient subcluster frequencies then plot by clinical group
plt_clinical <- Meta %>% 
  mutate(ReCluster = factor(annotated_SCs)) %>%           # convert subcluster labels to factor
  mutate(Groups = factor(Patient_groups, levels = clinical_groups)) %>%  # order clinical groups
  filter(ReCluster %in% subset_to_be_plotted) %>%         # keep only the selected CD4 subclusters
  group_by(Groups, Fig_ids, ReCluster) %>%                # group by clinical group, patient, and subcluster
  summarise(n = n()) %>% #, Set = first(Set)
  mutate(freq = n / sum(n) *100) %>%                      # compute each subcluster's % within the patient's CD4 T cells
  ungroup() %>%
  as.data.frame() %>% 
  #filter(ReCluster %in% subset_to_be_plotted) %>% 
  ggplot(aes(x = Groups, y = freq, fill = Groups, group = Groups)) +
  geom_boxplot(outlier.shape = NA) +                      # boxplot without outlier points (shown via jitter instead)
  geom_jitter(size = 0.2) +                               # overlay individual patient data points
  theme_bw()  +  #THEME +
  #ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +  # add significance brackets (t-test)
  #ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = F, vjust = 0.5) + 
  theme(legend.position = "none",                         # remove legend (groups are shown on x-axis)
        strip.text = element_text(size = 14, face='bold')) +
  facet_wrap(.~ReCluster, scales = "free_y", nrow = 1) +  # one panel per subcluster, independent y-axis scales
  
  scale_fill_manual(values=col_pGroups) + #**             # apply the predefined group colors
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle =0),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) + #    ylab('% PBMC') + xlab('Age groups')
  ylab('% in CD4 T cells') + xlab('Clinical groups')

print(plt_clinical)

