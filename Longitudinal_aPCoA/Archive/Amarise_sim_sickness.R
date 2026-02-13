load("mom_270.Rdata")
otu_original <- t(otu)

#libsize <- rowSums(otu_original)
p <- ncol(otu_original)
n = 100; m = 4

set.seed(30)
cond_taxa1 <- sample(1:p, 20)
cond_taxa2 <- sample(setdiff(1:p, cond_taxa1), 20)

set.seed(240)
selected_samples <- sample.int(nrow(otu_original), n)
cond <- rbinom(n, 1, 0.5)
#sick_time <- matrix(rbinom(n*m, 1, 0.5), m, n)
sick_time_options <- expand.grid(rep(list(0:1),m))
# sick_time_options <- sick_time_options[-which(apply(sick_time_options, 1, var) == 0),]
sick_time <- t(sick_time_options[sample(nrow(sick_time_options), n, replace = TRUE),])

# perturb independent draws of taxa counts for each sample
otu_tmp <- matrix(nrow = n*m, ncol = p+2)
otu_tmp <- as.data.frame(otu_tmp)
names(otu_tmp) <- c('patient.id','time',paste0('taxa',1:p))
otu_tmp$patient.id <- rep(1:n, each = m)
otu_tmp$time <- 1:m

# change condition and sick taxa
cond_effect <- c(0.25, 3, 2)
sick_effect <- c(24)
for (i in 1:n)
{
  otu_tmp[(i-1)*m+1, 3:(p+2)] <- bayesm::rdirichlet(otu_original[selected_samples[i],] + 0.5) *
    sum(otu_original[selected_samples[i],])
  if (sick_time[1,i])
  {
    otu_tmp[(i-1)*m+1, id + 2] <- otu_tmp[(i-1)*m+1, id + 2] * sick_effect[1]
  }
}

for (i in 1:n)
{
  for (j in 2:m)
  {
    # perturb previous time points' measurement to create current time point measurements
    otu_tmp[((i-1)*m + j), 3:(p+2)] <- bayesm::rdirichlet(as.numeric(otu_tmp[((i-1)*m + j - 1), 3:(p+2)]) + 0.5) *
      sum(otu_tmp[((i-1)*m + j - 1), 3:(p+2)])
    otu_tmp[((i-1)*m + j), cond_taxa1+2] <- otu_tmp[((i-1)*m + j), cond_taxa1+2] * cond_effect[1]
    
    if (cond[i])
    {
      otu_tmp[((i-1)*m + j), cond_taxa2+2] <- otu_tmp[((i-1)*m + j), cond_taxa2+2] * cond_effect[2]
    }
    
    if (sick_time[j, i])
    {
      # change id taxa
      otu_tmp[((i-1)*m + j), id + 2] = otu_tmp[((i-1)*m + j), id + 2] * sick_effect[1]
    }
  }
}

otu_tmp_unrounded <- otu_tmp
# round the otu table
otu_tmp <- round(otu_tmp)

# Aitchison
rel_tmp <- otu_tmp[,3:ncol(otu_tmp)]/rowSums(otu_tmp[,3:ncol(otu_tmp)])
rel_tmp_clr <- clr_lite(otu_tmp[,3:ncol(otu_tmp)], samples_are = "rows", method = "logunif", replicates = 100)

K <- D2K(as.matrix(dist(rel_tmp_clr)))
K_eigen <- eigen(K, symmetric = T)
mK <- which(K_eigen$values > 1e-9)
PCs <- K_eigen$vectors[ ,mK] %*% diag(sqrt(K_eigen$values[mK]))

sick <- as.vector(sick_time)

example_data <- list(metadata = cbind(otu_tmp[,1:2],
                                      batch = sick,
                                      treatment = rep(cond, each = m)) %>%
                       rename(subjectid = `patient.id`) %>%
                       # SL 1.18:
                       mutate(time = factor(time,
                                            levels = sort(unique(time)),
                                            labels = paste0("t", sort(unique(time)))
                       )),
                     otu_counts = otu_tmp[,3:ncol(otu_tmp)])

# add sick as a fixed effect
mK90 <- min(which(cumsum(K_eigen$values[mK]/sum(K_eigen$values[mK])) > 0.9))
PCs <- K_eigen$vectors[ ,1:mK90] %*% diag(sqrt(K_eigen$values[1:mK90]))
colnames(PCs) <- paste0("PC", seq_len(ncol(PCs)))


# Linear mixed model
m1 <- lmm_long_apcoa(PCs, 
                     otu_tmp, 
                     sick, 
                     n, 
                     m)
resid_PCs <- lapply(m1, function(x) process_resids(x, otu_tmp, cond, sick))
p_final1 <- plot_apcoa(resid_PCs[[3]], "LMM")
p_final1



# Quantile
m2 <- quantile_apcoa(PCs = PCs, 
                     metadata = example_data$metadata, 
                     to_remove = c("batch","time"),
                     to_remain = c("treatment"),
                     taus = seq(0.01, 0.99, length.out = 99),
                     tauw = rep(1/99, 99),
                     subject_id = "subjectid"
)
final_PCs <- process_resids(m2, otu_tmp, cond, sick)
p_final2 <- plot_apcoa(final_PCs,"Quantile")
p_final2

