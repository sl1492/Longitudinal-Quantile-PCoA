## Simulation Setup
simulate_data <- function(rdata_path = "mom_270.Rdata",
                          n = 100, m = 4,
                          cond_effect = 0.25, batch_effect = 8,
                          seed_init = 1, seed_select = 7711) {
  
  load(rdata_path)
  otu_original <- t(otu)
  p <- ncol(otu_original)
  n = 100 
  m = 4

  # change condition taxa
  set.seed(1)
  cond_taxa1 <- sample(setdiff(1:p, id), 20) # introduce changes in the community profiles shared by all subjects
  cond_taxa2 <- sample(setdiff(1:p, cond_taxa1), 20) # Treatment group, could have both sex and treatment effects

  set.seed(1)
  selected_samples <- sample.int(nrow(otu_original), n)

  cond <- rbinom(n, 1, 0.5) # treatment/control
  batchid <- rbinom(n, 1, 0.5) # sex

  # perturb independent draws of taxa counts for each sample
  otu_tmp <- matrix(nrow = n*m, ncol = p+2)
  otu_tmp <- as.data.frame(otu_tmp)
  names(otu_tmp) <- c('patient.id','time',paste0('taxa',1:p))
  otu_tmp$patient.id <- rep(1:n, each = m)
  otu_tmp$time <- 1:m

  # Initialize taxa counts
  for (i in 1:n) {
    otu_tmp[(i-1)*m + 1, 3:(p+2)] <- bayesm::rdirichlet(otu_original[selected_samples[i],] + 0.5) *
      sum(otu_original[selected_samples[i],])
  }
  
  cond_effect <- c(0.25) # introduce changes in the community profiles shared by all subjects
  batch_effect <- c(8) # obscuring sex effect
  sub_fc <- ifelse(cond == 1, 6 * rlnorm(n, meanlog = log(5), sdlog = 1.2), 1) # SL 10.17

  for (i in 1:n)
  {
    for (j in 2:m)
    {
      # perturb previous time points' measurement to create current time point measurements
      otu_tmp[((i-1)*m + j), 3:(p+2)] <- bayesm::rdirichlet(as.numeric(otu_tmp[((i-1)*m + j - 1), 3:(p+2)])) *
        sum(otu_tmp[((i-1)*m + j - 1), 3:(p+2)])
      otu_tmp[((i-1)*m + j), (cond_taxa1 + 2)] <- otu_tmp[((i-1)*m + j), (cond_taxa1 + 2)] * cond_effect[1]
      
      if (cond[i]) # for treated subjects only
      {
        otu_tmp[((i-1)*m + j), (cond_taxa2 + 2)] <- otu_tmp[((i-1)*m + j), (cond_taxa2 + 2)] * sub_fc[i]
      }
      if (batchid[i]) # for female only
      {
        otu_tmp[((i-1)*m + j), (id + 2)] = otu_tmp[((i-1)*m + j), (id + 2)] * batch_effect[1]
      }
    }
  }
  otu_tmp_unrounded <- otu_tmp
  otu_tmp <- round(otu_tmp)
  
  example_data <- list(metadata = cbind(otu_tmp[,1:2],
                                        batch = rep(batchid, each = m),
                                        treatment = rep(cond, each = m)) %>%
                         rename(subjectid=`patient.id`),
                       otu_counts = otu_tmp[,3:ncol(otu_tmp)])
  
  list(
    otu_tmp = otu_tmp,
    otu_tmp_unrounded = otu_tmp_unrounded,
    example_data = example_data,
    cond = cond,
    batchid = batchid,
    n = n, m = m, p = p
  )
}
