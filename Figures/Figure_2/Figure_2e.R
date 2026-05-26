# Load required libraries
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)


# Load metadata containing PBMC cell information (n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# Define color palette for the four CD14 monocyte subclusters
cols <- c('CD14_mono_SC0'='tomato','CD14_mono_SC1'='paleturquoise','CD14_mono_SC2'='cornflowerblue','CD14_mono_SC3'='mediumseagreen')

# Specify which monocyte subclusters to include in the plot
subset_to_be_plotted <- c( 'CD14_mono_SC0','CD14_mono_SC1','CD14_mono_SC2', 'CD14_mono_SC3')

# Define individual groups: 14 healthy controls (pHC) and 26 COVID-19 patients grouped by severity (G1, G2, G3)
pHC <-  paste0("pHC", seq(1:14))
pG1 <- c("pCoV1", "pCoV4", "pCoV10", "pCoV13", "pCoV14", "pCoV15", "pCoV16", "pCoV18", "pCoV19", "pCov25")
pG2 <-  c("pCoV3",  "pCoV5",  "pCoV6",  "pCoV8", "pCoV9", "pCoV12", "pCoV20", "pCov22" ,"pCov23", "pCov24", "pCov26")
pG3 <- c( "pCoV2","pCoV7", "pCoV11", "pCoV17","pCov21")
ordered_names <- c(pHC, pG1, pG2, pG3)

# Create stacked bar plot showing frequency of each monocyte subcluster per individual
BP <- Meta %>% 
  
  # Convert Groups to ordered factor for consistent visualization
  mutate(Patient_groups = factor(Groups, levels = c("pHC", "G1", "G2", "G3"))) %>%
  # Convert annotated_SCs to factor for grouping
  mutate(ReCluster = factor(annotated_SCs)) %>%
  #mutate(ReCluster = factor(Names, levels = order_names)) %>%
  # Filter to include only the specified monocyte subclusters
  filter(ReCluster %in% subset_to_be_plotted) %>%
  
  # Count cells per individual per subcluster per group
  group_by(Groups, Fig_ids, ReCluster) %>%
  #filter(Groups %in% c("HO_M",'HO_F')) %>% 
  summarise(n = n()) %>% #, Set = first(Set)
  # Calculate frequency (percentage) within each individual
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>%
  
  # Generate stacked bar plot
  ggplot(aes(x = Fig_ids, y = freq, fill = ReCluster, group = Groups)) +
  geom_bar(stat = "identity") +
  # Apply custom color palette
  scale_fill_manual(values=cols) + #***
  # Order x-axis by individual groups (healthy controls first, then patient groups)
  scale_x_discrete(limits=ordered_names) + #labels= labels
  # Customize plot aesthetics
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle = 90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") + #    ylab('% PBMC') + xlab('Age groups')
  
  # Add axis labels
  ylab('% in CD14 mo') + xlab('Individuals (n=40)')

# Display the plot
print(BP)
