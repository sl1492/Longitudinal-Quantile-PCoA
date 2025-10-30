library(bindata)
library(GENESIS)
library(dplyr)
library(ggplot2)
library(Matrix)
library(nlme)
library(lme4)
library(Tjazi)
library(cowplot)
library(rqpd)
library(quantreg)
library(zCompositions)
library(compositions)

## Helper functions

# SL 9.26: from MiRKAT
D2K <- function(D) {
  
  n <- nrow(D)
  centerM <- diag(n) - 1/n*(rep(1,n)%*%t(rep(1,n)))
  K <- -0.5 * centerM %*% (D * D) %*% centerM # Gower Centering
  eK <- eigen(K, symmetric = TRUE)
  K <- eK$vector %*% diag(pmax(0, eK$values)) %*% t(eK$vector) # SL 10.5: pmax(0,..) ensure PSD
  
  return(K)
}

# Function to produce K* (last step, adapted from Amarise's code)
process_resids <- function(resids, 
                           otu_tmp, 
                           cond, 
                           batchid) {
  
  eigen <- eigen(tcrossprod(resids), symmetric = TRUE)
  mK <- which(eigen$values > 1e-9) # keep positive values
  PCs <- eigen$vectors[,mK] %*% diag(sqrt(eigen$values[mK]))/1
  PCs <- as.data.frame(cbind(patient.id = otu_tmp$patient.id,
                             time = otu_tmp$time,
                             treatment = rep(cond, each = max(otu_tmp$time)),
                             batch = rep(batchid, each = max(otu_tmp$time)),
                             PCs))
  names(PCs)[5:ncol(PCs)] <- c(paste0('PC',1:(ncol(PCs)-4)))
  
  return(list(PCs = PCs,
              PC1_perc = eigen$values[1]/sum(eigen$values[mK]),
              PC2_perc = eigen$values[2]/sum(eigen$values[mK])))
}

# The ploting function
plot_apcoa <- function(pc_list, title, pc1_name = "PC1", pc2_name = "PC2",
                     pc1_perc = NULL, pc2_perc = NULL) {
  
  annotations <- data.frame(
    time = 1:4,
    label = c("t[0]","t[1]","t[2]","t[3]"),
    x = -Inf, y = Inf
  )
  
  pc1_lab <- if (is.null(pc1_perc)) pc1_name else paste0(pc1_name, " (", round(100*pc1_perc, 2), "%)")
  pc2_lab <- if (is.null(pc2_perc)) pc2_name else paste0(pc2_name, " (", round(100*pc2_perc, 2), "%)")
  
  pc_list$PCs %>%
    mutate(Sex = ifelse(batch == 1, "Female", "Male"),
           Arm = factor(ifelse(treatment == 0, "Control", "Treatment"),
                        levels = c("Treatment","Control"))) %>%
    ggplot(aes(PC1, PC2, color = Arm, shape = Sex, group = Sex)) +
    geom_point(size = 2) +
    facet_wrap(~time, nrow = 2, scales = "free",
               labeller = label_bquote(t[.(time-1)])) +
    theme(strip.text = element_blank()) +
    xlab(pc1_lab) + ylab(pc2_lab) +
    theme(plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18)) +
    ggtitle(title) +
    stat_ellipse(aes(linetype = Sex)) +
    geom_text(data = annotations, aes(x = x, y = y, label = label),
              hjust = -0.2, vjust = 1.5, parse = TRUE, inherit.aes = FALSE)
}
