
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyr)


 #' A function for principal variance component analysis
   #'
   #' The function is written based on the 'pvcaBatchAssess' function of the PVCA R package
   #' and slightly changed to make it more efficient and flexible for sequencing read counts data.
   #' (http://watson.nci.nih.gov/bioc_mirror/packages/release/bioc/manuals/pvca/man/pvca.pdf)
   #'
   #' @param counts The Normalized(e.g. TMM)/ log-transformed reads count matrix from sequencing data (row:gene/feature, col:sample)
   #' @param meta  The Meta data matrix containing predictor variables (row:sample, col:predictor)
   #' @param threshold The proportion of the variation in read counts explained by top k PCs. This value determines the number of top PCs to be used in pvca.
   #' @param inter TRUE/FALSE - include/do not include pairwise interactions of predictors
   #'
   #' @return std.prop.val The vector of proportions of variation explained by each predictor.
   #'
   #' @export
   #'
   
   PVCA <- function(counts, meta, threshold, inter){
     
     counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
     cor.counts <- cor(counts.center)
     dim(cor.counts)
     eigen.counts <- eigen(cor.counts)
     eigen.mat <- eigen.counts$vectors
     eigen.val <- eigen.counts$values
     n.eigen <- length(eigen.val)
     eigen.val.sum <- sum(eigen.val)
     percents.pcs <- eigen.val/eigen.val.sum
     meta <- as.data.frame(meta)
     
     all <- 0
     npc.in <- 0
     for(i in 1:n.eigen){
       all <- all + percents.pcs[i]
       npc.in <- npc.in + 1
       if(all > threshold){break}
     }
     if (npc.in < 3) {npc <- 3}
     
     pred.list <- colnames(meta)
     meta <- droplevels(meta)
     
     n.preds <- ncol(meta) + 1
     if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}
     
     ran.pred.list <- c()
     for(i in 1:ncol(meta)){
       ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
     }
     ##interactions
     if(inter){
       for(i in 1:(ncol(meta)-1)){
         for(j in (i+1):ncol(meta)){
           ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
           pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
         }
       }
     }
     formula <- paste(ran.pred.list, collapse = " + ")
     formula <- paste("pc", formula, sep=" ~ ")
     ran.var.mat <- NULL
     for(i in 1:npc.in){
       dat <- cbind(eigen.mat[,i],meta)
       colnames(dat) <- c("pc",colnames(meta))
       Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit)
       var.vec <- unlist(VarCorr(Rm1ML))
       ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
     }
     ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
     wgt.vec <- eigen.val/eigen.val.sum
     prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
     std.prop.var <- prop.var/sum(prop.var)
     std.prop.var
   }
   
   # load sample information   
   meta <- read.table('pCov40_SampleInfo_03112025.txt', header = T, row.names = 1) 
   head(meta)
   meta <- meta[,-c(1,2,3)]
   
   #Cell compositions (sample information)
   y=prop.table(x=table(Meta$IDs, as.character(Meta$SCs)), margin=2)  
   head(y)
   counts = t(y)
   head(counts)
  
   
   # Reorder rows in metadata to align with column names in counts
   meta <- meta[match(colnames(counts), rownames(meta)), ]
   dim(meta)
   head(meta)
   
   #set threshold
   threshold = 0.2
   
   #set inter 
   inter = FALSE 
   
  # Run PVCA
   results <- PVCA(counts, meta, threshold, inter)
   print(results)
   
  # plot 
   PlotPVCA <- function(pvca.res, title){
     plot.dat <- data.frame(eff=names(pvca.res), prop=pvca.res)
     p <- ggplot2::ggplot(plot.dat, aes(x=eff, y=prop))
     p <- p + ggplot2::ggtitle(title)
     p <- p + ggplot2::geom_bar(stat="identity", fill="#bcbddc", colour="#bcbddc")
     p <- p + ggplot2::geom_text(aes(label=round(prop,3), y=prop+0.04), size=4)
     p <- p + ggplot2::scale_x_discrete(limits=names(pvca.res))
     p <- p + ggplot2::scale_y_continuous(limits = c(0,1))
     p <- p + ggplot2::labs(x= "Effects", y= "Weighted average proportion variance")
     p <- p + ggplot2::theme_bw()
     p <- p + ggplot2::theme(plot.background = element_blank() ,panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank() ,panel.border = element_blank(), panel.background = element_blank())
     p <- p + ggplot2::theme(axis.line = element_line(color = 'black'))
     p <- p + ggplot2::theme(axis.title.x = element_text(size = 15, vjust=-0.5))
     p <- p + ggplot2::theme(axis.title.y = element_text(size = 15, vjust= 1.0))
     p <- p + ggplot2::theme(axis.text = element_text(size = 12))
     p <- p + ggplot2::theme(axis.text.x = element_text(angle = 90, vjust= 0.5, hjust=1))
     p
   }
   
   plot <- PlotPVCA(results, '')
   plot
   
   
   ggsave("./PVCA_analysis_pCoV_12022024.pdf", plot,
          width=7, height=3.8,  units="in", scale=1)
   
