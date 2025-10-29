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

## Construct Aitchison Kernel and extract PCs

# Construct D
rel_tmp <- otu_tmp[,3:ncol(otu_tmp)]/rowSums(otu_tmp[,3:ncol(otu_tmp)]) # may give NA values if denominator is 0

# SL 9.26: Compute the corresponding CLR (centered log-ratio) transformed abundance (Amarise used clr_lite)
imputed = cmultRepl(rel_tmp, suppress.print = TRUE) # impute zero
rel_tmp_clr = as.matrix(clr(imputed))
rownames(rel_tmp_clr) <- rownames(imputed)

# SL 10.17: since cmultRepl() can drop NA rows, the following block is added
shared <- intersect(rownames(rel_tmp_clr), rownames(example_data$metadata))
rel_tmp_clr <- rel_tmp_clr[shared, ]
metadata <- example_data$metadata[shared, ]

# Construct K
K <- D2K(as.matrix(dist(rel_tmp_clr)))
K_eigen <- eigen(K, symmetric = T)

# SL 9.26: Adapted Amarise's code
# set of positive eigenvalues. For a PSD kernel, the count of positive eigenvalues equals rank(K)
mK <- which(K_eigen$values > 1e-9) 

# Retain the top ℓ ≤ rank(K) kernel PCs that explain a large proportion of the variability in K
mK90 <- min(which(cumsum(K_eigen$values[mK]/sum(K_eigen$values[mK])) > 0.9))
PCs <- K_eigen$vectors[ ,1:mK90] %*% diag(sqrt(K_eigen$values[1:mK90]))
