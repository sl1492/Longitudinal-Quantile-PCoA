source("00_helpers.R")
source("01_simulation.R")
source("02_kernel_pca.R")
source("03_linear_aPCoA.R")
source("04_quantile_aPCoA.R")

# simulate data
sim <- simulate_data("mom_270.Rdata", n = 100, m = 4)

# build Aitchison kernel & PCs
kpca <- build_kernel_pcs(sim$otu_tmp, sim$example_data)

# LMM
m1 <- lmm_long_apcoa(kpca$PCs, 
                     sim$otu_tmp, 
                     sim$batchid, 
                     sim$cond, 
                     sim$n, 
                     sim$m)
p1 <- plot_apcoa(m1[[3]], "LMM")

# Quantile
m2 <- quantile_apcoa(PCs = kpca$PCs, 
                     metadata = kpca$metadata, 
                     covariates = "all", 
                     lambda = 1.5, q = 99, 
                     formula = "batch + time",
                     otu_tmp = sim$otu_tmp, 
                     batchid = sim$batchid, 
                     cond = sim$cond
                     )
p2 <- plot_apcoa(m2,"Quantile")

# plots

p1; p2  # show in Viewer
