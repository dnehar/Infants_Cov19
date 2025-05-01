




# colors  
col_SCs <- c('CD16_mono_SC0'='paleturquoise',
             'CD16_mono_SC1'='tomato',
             'CD16_mono_SC2'='cornflowerblue',
             'CD16_mono_SC3'='mediumseagreen')


# B cell subsets 
subset_to_be_plotted <- c("CD16_mono_SC0", paste0("CD16_mono_SC",seq(1:3)))

# order names 
pHC <-  paste0("pHC", seq(1:14))
pG1 <- c("pCoV1", "pCoV4", "pCoV10", "pCoV13", "pCoV14", "pCoV15", "pCoV16", "pCoV18", "pCoV19", "pCov25")
pG2 <-  c("pCoV3",  "pCoV5",  "pCoV6",  "pCoV8", "pCoV9", "pCoV12", "pCoV20", "pCov22" ,"pCov23", "pCov24", "pCov26")
pG3 <- c( "pCoV2","pCoV7", "pCoV11", "pCoV17","pCov21")
aHC <-  c(paste0("aHC", seq(1:3)), paste0("aHC", seq(7, 14, by=1))) 
aG1 <- c("aCoV8","aCoV9","aCoV20")
aG2 <- c("aCoV3", "aCoV21", "aCoV25", "aCoV28"  )
aG3 <- c("aCoV1", "aCoV2" ,"aCoV4",  "aCoV5",  "aCoV6",   "aCoV7", "aCoV10", "aCoV11", "aCoV12", "aCoV13",
         "aCoV14", "aCoV15", "aCoV16", "aCoV17", "aCoV18", "aCoV19",  "aCoV22", "aCoV23", "aCoV24", "aCoV26", "aCoV27") 

ordered_names <- c(pHC, pG1, pG2, pG3, aHC, aG1, aG2, aG3)
length(ordered_names)


BP <- Meta %>% 
  
  mutate(Patient_groups = factor(Groups, levels = clinical_groups)) %>%
  mutate(ReCluster = factor(SCs)) %>%
  #mutate(ReCluster = factor(Names, levels = order_names)) %>%
  filter(ReCluster %in% subset_to_be_plotted) %>%
  
  group_by(Groups, Fig_ids, ReCluster) %>%
  #filter(Groups %in% c("HO_M",'HO_F')) %>% 
  summarise(n = n()) %>% #, Set = first(Set)
  mutate(freq = n / sum(n) *100) %>%
  ungroup() %>%
  as.data.frame() %>%
  ggplot(aes(x = Fig_ids, y = freq, fill = ReCluster, group = Groups)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=col_SCs) + #***
  scale_x_discrete(limits=ordered_names) + #labels= labels
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_text(size=16, angle = 90),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") + #    ylab('% PBMC') + xlab('Age groups')
  
  ylab('% in CD16 mo') + xlab('Individuals (n=40)')

print(BP)
