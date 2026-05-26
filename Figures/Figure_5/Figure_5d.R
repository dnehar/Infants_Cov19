# Load required libraries for data manipulation and visualization
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)


# Load metadata file containing B cell subcluster information (PBMCs, n=203,402 cells)   
Meta <- read.csv('Meta_pCoV40_03112025_small.csv')

# Define colors for each B cell subcluster for consistent visualization
cols <- c('B_cells_SC0'='tomato','B_cells_SC1'='paleturquoise','B_cells_SC3'='mediumseagreen','B_cells_SC2'='cornflowerblue','B_cells_SC4'='mediumpurple')

# Specify which B cell subclusters to include in the plot
subset_to_be_plotted <- c('B_cells_SC0','B_cells_SC1','B_cells_SC2','B_cells_SC3','B_cells_SC4')

# Define donor groups for ordering samples on x-axis
# pHC: healthy controls (n=14), pG1-pG3: COVID-19 patient groups stratified by some criteria
pHC <-  paste0("pHC", seq(1:14))
pG1 <- c("pCoV1", "pCoV4", "pCoV10", "pCoV13", "pCoV14", "pCoV15", "pCoV16", "pCoV18", "pCoV19", "pCov25")
pG2 <-  c("pCoV3",  "pCoV5",  "pCoV6",  "pCoV8", "pCoV9", "pCoV12", "pCoV20", "pCov22" ,"pCov23", "pCov24", "pCov26")
pG3 <- c("pCoV2","pCoV7", "pCoV11", "pCoV17","pCov21")
ordered_names <- c(pHC, pG1, pG2, pG3)

# Create and print stacked bar plot
BP <- Meta %>% 
  # Convert Groups to factor with specific order for consistent grouping
  mutate(Patient_groups = factor(Groups, levels = c("pHC", "G1", "G2", "G3"))) %>%
  # Convert annotated subclusters to factor
  mutate(ReCluster = factor(annotated_SCs)) %>%
  # Keep only the B cell subclusters specified above
  filter(ReCluster %in% subset_to_be_plotted) %>%
  # Count cells per group, individual, and subcluster
  group_by(Groups, Fig_ids, ReCluster) %>%
  summarise(n = n()) %>%
  # Calculate frequency (percentage) of each subcluster within each individual
  mutate(freq = n / sum(n) *100) %>%
  # Ungroup and convert to dataframe for plotting
  ungroup() %>%
  as.data.frame() %>%
  # Create stacked bar plot with individuals on x-axis and percentage on y-axis
  ggplot(aes(x = Fig_ids, y = freq, fill = ReCluster, group = Groups)) +
  geom_bar(stat = "identity") +
  # Apply custom colors for each B cell subcluster
  scale_fill_manual(values=cols) +
  # Order individuals on x-axis according to predefined order (HC first, then patient groups)
  scale_x_discrete(limits=ordered_names) +
  # Customize theme: increase font sizes and rotate x-axis labels for readability
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle = 90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") +
  # Label axes
  ylab('% in B cells') + xlab('Individuals (n=40)')

# Display the final plot
print(BP)
