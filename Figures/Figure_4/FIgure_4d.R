# Load required libraries for data manipulation and plotting
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Load metadata file containing PBMC single-cell data (203,402 cells)
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# Define color mapping for each CD8 T cell subcluster
cols <- c('CD8_NAIVE'='tomato','CD8_ISGhi'='paleturquoise','CD8_GzK'='mediumseagreen','CD8_TEMRA'='cornflowerblue','CD8_Prolif'='mediumpurple')

# Specify which CD8 T cell subclusters to include in the plot
subset_to_be_plotted <- c( 'CD8_NAIVE','CD8_ISGhi','CD8_GzK','CD8_TEMRA','CD8_Prolif')

# Define sample grouping: healthy controls (pHC), and three patient groups (pG1, pG2, pG3)
pHC <-  paste0("pHC", seq(1:14))
pG1 <- c("pCoV1", "pCoV4", "pCoV10", "pCoV13", "pCoV14", "pCoV15", "pCoV16", "pCoV18", "pCoV19", "pCov25")
pG2 <-  c("pCoV3",  "pCoV5",  "pCoV6",  "pCoV8", "pCoV9", "pCoV12", "pCoV20", "pCov22" ,"pCov23", "pCov24", "pCov26")
pG3 <- c( "pCoV2","pCoV7", "pCoV11", "pCoV17","pCov21")
# Combine all samples in desired order
ordered_names <- c(pHC, pG1, pG2, pG3)

# Create stacked bar plot
BP <- Meta %>% 
  # Convert Groups to ordered factor for consistent legend ordering
  mutate(Patient_groups = factor(Groups, levels = c("pHC", "G1", "G2", "G3"))) %>%
  # Convert cell cluster annotations to factor
  mutate(ReCluster = factor(annotated_SCs)) %>%
  # Filter to include only selected CD8 T cell subclusters
  filter(ReCluster %in% subset_to_be_plotted) %>%
  
  # Group data by patient group, individual (Fig_ids), and cell cluster
  group_by(Groups, Fig_ids, ReCluster) %>%
  # Count cells in each group
  summarise(n = n()) %>%
  # Calculate frequency (percentage) within each individual
  mutate(freq = n / sum(n) *100) %>%
  # Remove grouping
  ungroup() %>%
  as.data.frame() %>%
  
  # Create ggplot visualization with x=individuals, y=percentage, stacked by cluster
  ggplot(aes(x = Fig_ids, y = freq, fill = ReCluster, group = Groups)) +
  # Use stacked bar chart
  geom_bar(stat = "identity") +
  # Apply predefined colors to subclusters
  scale_fill_manual(values=cols) +
  # Order x-axis according to sample groups (HC first, then G1, G2, G3)
  scale_x_discrete(limits=ordered_names) +
  # Set theme: increase font sizes, rotate x-axis labels, remove legend
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle = 90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") +
  # Set axis labels
  ylab('% in CD8 T cells') + xlab('Individuals (n=40)')

# Display the plot
print(BP)
