## Construct Aitchison Kernel and extract PCs
build_kernel_pcs <- function(otu_tmp, example_data) {
  
  # Construct D
  rel_tmp <- otu_tmp[,3:ncol(otu_tmp)]/rowSums(otu_tmp[,3:ncol(otu_tmp)]) # may give NA values if denominator is 0
  
  # SL 10.20: Filter out the low abundant genera from the relative abundance data
  thres_abundance = 0.01 / 100
  mean_abundance <- apply(rel_tmp, 2, function(x) mean(x, na.rm = TRUE))
  rel_tmp <- rel_tmp[, mean_abundance > thres_abundance]
  
  # SL 9.26: Compute the corresponding CLR (centered log-ratio) transformed abundance 
  imputed = cmultRepl(rel_tmp, suppress.print = TRUE) # impute zero
  rel_tmp_clr = as.matrix(clr(imputed))
  rownames(rel_tmp_clr) <- rownames(imputed)
  
  # Construct K
  K <- D2K(as.matrix(dist(rel_tmp_clr)))
  K_eigen <- eigen(K, symmetric = T)
  
  # SL 9.26: Adapted Amarise's code
  # set of positive eigenvalues. For a PSD kernel, the count of positive eigenvalues equals rank(K)
  mK <- which(K_eigen$values > 1e-9) 
  
  # Retain the top ℓ ≤ rank(K) kernel PCs that explain a large proportion of the variability in K
  mK90 <- min(which(cumsum(K_eigen$values[mK]/sum(K_eigen$values[mK])) > 0.9))
  PCs <- K_eigen$vectors[ ,1:mK90] %*% diag(sqrt(K_eigen$values[1:mK90]))
  colnames(PCs) <- paste0("PC", seq_len(ncol(PCs)))
  
  list(PCs = PCs,
       K_eigen = K_eigen, 
       mK90 = mK90)
}

