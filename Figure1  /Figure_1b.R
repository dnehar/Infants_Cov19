
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)


VL <- read.csv('pCoV_2022/pCoV40_ViralLoads.csv')
head(VL)  
table(VL$Groups)

p_VL <- ggplot(VL, aes(x=Groups, y=VL, fill=Groups)) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ #+ #binwidth = 1
  ylab("SARS-CoV-2 PCR CT values") +
  xlab("Groups") +
  scale_y_continuous(limits = c(0,45)) +
  #geom_hline(aes(yintercept=30), color="black",linetype="dashed")+
  scale_fill_manual(values=col_pGroups) + 
  theme(panel.background = element_rect(fill = "white", colour = NA),
        axis.title.x = element_text(face="bold", size=0),
        axis.title.y = element_text(face="bold", size=14),
        axis.text.y=element_text(face="bold",size=12), 
        axis.text.x=element_text(face="bold",size=12), 
        # plot.title = element_text(hjust = 0.5,face='bold',size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, color = 'black', size=1),
        legend.position = "none")
#ggtitle('Patient groups', )
p_VL 
