# Load required libraries
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Define color palette for the four NK cell subsets
cols <- c('NK_SC0'='tomato','NK_SC1'='paleturquoise','NK_SC3'='mediumseagreen','NK_SC2'='cornflowerblue')

# Specify which NK cell subsets to include in the plot
subset_to_be_plotted <- c( 'NK_SC0','NK_SC1','NK_SC2','NK_SC3')

# Create vectors of patient/sample IDs for each group
pHC <- paste0("pHC", seq(1:14))  # 14 healthy controls
pG1 <- c("pCoV1", "pCoV4", "pCoV10", "pCoV13", "pCoV14", "pCoV15", "pCoV16", "pCoV18", "pCoV19", "pCov25")  # 10 samples in group 1
pG2 <- c("pCoV3",  "pCoV5",  "pCoV6",  "pCoV8", "pCoV9", "pCoV12", "pCoV20", "pCov22" ,"pCov23", "pCov24", "pCov26")  # 11 samples in group 2
pG3 <- c( "pCoV2","pCoV7", "pCoV11", "pCoV17","pCov21")  # 5 samples in group 3

# Combine all IDs in order for x-axis plotting
ordered_names <- c(pHC, pG1, pG2, pG3)

# Create and plot the bar chart
BP <- Meta %>% 
  # Convert Groups to ordered factor (pHC, G1, G2, G3)
  mutate(Patient_groups = factor(Groups, levels = c("pHC", "G1", "G2", "G3"))) %>%
  # Convert annotated_SCs to factor for NK cell subset grouping
  mutate(ReCluster = factor(annotated_SCs)) %>%
  # Keep only the specified NK cell subsets
  filter(ReCluster %in% subset_to_be_plotted) %>%
  
  # Group by patient group, individual ID, and NK subset
  group_by(Groups, Fig_ids, ReCluster) %>%
  # Count number of cells in each group
  summarise(n = n()) %>%
  # Calculate percentage frequency within each individual
  mutate(freq = n / sum(n) *100) %>%
  # Ungroup to prepare for plotting
  ungroup() %>%
  as.data.frame() %>%
  # Generate stacked bar plot
  ggplot(aes(x = Fig_ids, y = freq, fill = ReCluster, group = Groups)) +
  geom_bar(stat = "identity") +
  # Apply custom color palette
  scale_fill_manual(values=cols) +
  # Order x-axis by predefined patient order
  scale_x_discrete(limits=ordered_names) +
  # Customize plot appearance (text size, angle, bold labels)
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle = 90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") +
  # Set axis labels
  ylab('% in NK cells') + xlab('Individuals (n=40)')

# Display the plot
BP 
