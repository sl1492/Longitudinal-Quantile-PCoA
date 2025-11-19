## Simulation Setup
simulate_data <- function(rdata_path = "mom_270.Rdata",
                          n = 100, m = 4,
                          cond_effect = 4, batch_effect = 8,
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

  # SL 11.14: split into decreasing and increasing for balanced effect
  id_d <- 1:10    # first 10 decrease
  id_i <- 11:20   # last 10 increase
  
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
  
  cond_effect <- c(4) # introduce changes in the community profiles shared by all subjects
  batch_effect <- c(8) # obscuring sex effect
  treat_effect <- c(3) # treatment effect
  #sub_fc <- ifelse(cond == 1, 6 * rlnorm(n, meanlog = log(5), sdlog = 1.2), 1) # SL 10.17 treatment effect

  for (i in 1:n)
  {
    for (j in 2:m)
    {
      # SL 11.14
      prev_row <- (i-1)*m + (j-1)
      cur_row  <- (i-1)*m + j
      
      # perturb previous time points' measurement to create current time point measurements
      prev_counts <- as.numeric(otu_tmp[prev_row, 3:(p+2)])
      otu_tmp[cur_row, 3:(p+2)] <- bayesm::rdirichlet(prev_counts + 0.5) * sum(prev_counts)
      # otu_tmp[cur_row, (cond_taxa1 + 2)] <- otu_tmp[cur_row, (cond_taxa1 + 2)] * cond_effect[1]
      
      # SL 11.19
      s1_cond <- sum(otu_tmp[cur_row, cond_taxa1[id_d] + 2], na.rm = TRUE)
      s2_cond <- sum(otu_tmp[cur_row, cond_taxa1[id_i] + 2], na.rm = TRUE)

      if (s1_cond > 0 && s2_cond > 0) {
        d <- cond_effect[1]
        otu_tmp[cur_row, cond_taxa1[id_d] + 2] <- otu_tmp[cur_row, cond_taxa1[id_d] + 2] / d
        change_fold <- 1 + (s1_cond/s2_cond)*((d - 1)/d)
        otu_tmp[cur_row, cond_taxa1[id_i] + 2] <- otu_tmp[cur_row, cond_taxa1[id_i] + 2] * change_fold
      }
      
      if (cond[i]) # for treated subjects only
      {
        # SL 11.19
        s1_treat <- sum(otu_tmp[cur_row, cond_taxa2[id_d] + 2], na.rm = TRUE)
        s2_treat <- sum(otu_tmp[cur_row, cond_taxa2[id_i] + 2], na.rm = TRUE)

        if (s1_treat > 0 && s2_treat > 0) {
          d <- treat_effect[1] 
          otu_tmp[cur_row, cond_taxa2[id_d] + 2] <- otu_tmp[cur_row, cond_taxa2[id_d] + 2] / d
          change_fold <- 1 + (s1_treat/s2_treat)*((d - 1)/d)
          otu_tmp[cur_row, cond_taxa2[id_i] + 2] <- otu_tmp[cur_row, cond_taxa2[id_i] + 2] * change_fold
        }
      }
      
      if (batchid[i]) # for female only
      {
        # SL 11.19
        s1_batch <- sum(otu_tmp[cur_row, id[id_d] + 2], na.rm = TRUE)
        s2_batch <- sum(otu_tmp[cur_row, id[id_i] + 2], na.rm = TRUE)
        
        if (s1_batch > 0 && s2_batch > 0) {
          d <- batch_effect[1]
          otu_tmp[cur_row, id[id_d] + 2] <- otu_tmp[cur_row, id[id_d] + 2] / d
          change_fold <- 1 + (s1_batch/s2_batch)*((d - 1)/d)
          otu_tmp[cur_row, id[id_i] + 2] <- otu_tmp[cur_row, id[id_i] + 2] * change_fold       
        }
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
