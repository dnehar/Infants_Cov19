# Load required libraries for data manipulation and visualization
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)

# Load viral load data from CSV file
VL <- read.csv('pCoV_2022/pCoV40_ViralLoads.csv')
head(VL)  # Display first few rows of data
table(VL$Groups)  # Show count of samples in each group

# Create violin plot with individual data points overlaid
p_VL <- ggplot(VL, aes(x=Groups, y=VL, fill=Groups)) + 
  geom_violin() +  # Draw distribution curves for each group
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ # Overlay individual data points
  ylab("SARS-CoV-2 PCR CT values") +  # Label y-axis
  xlab("Groups") +  # Label x-axis
  scale_y_continuous(limits = c(0,45)) +  # Set y-axis range from 0 to 45
  #geom_hline(aes(yintercept=30), color="black",linetype="dashed")+  # Optional: add reference line at CT=30
  scale_fill_manual(values=col_pGroups) +  # Apply custom color scheme to groups
  theme(panel.background = element_rect(fill = "white", colour = NA),  # White background
        axis.title.x = element_text(face="bold", size=0),  # Hide x-axis title
        axis.title.y = element_text(face="bold", size=14),  # Bold y-axis title
        axis.text.y=element_text(face="bold",size=12),  # Bold y-axis tick labels
        axis.text.x=element_text(face="bold",size=12),  # Bold x-axis tick labels
        # plot.title = element_text(hjust = 0.5,face='bold',size=14),  # Optional: add title
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.border = element_rect(fill=NA, color = 'black', size=1),  # Add black border
        legend.position = "none")  # Hide legend
#ggtitle('Patient groups', )  # Optional: add title
p_VL  # Display the plot
