### Simulation Setup
simulate_data_sex <- function(rdata_path = "mom_270.Rdata",
                          n = 100, # subjects
                          m = 4, # time point
                          batch_effect = c(1.25), # obscuring sex effect
                          treat_effect = c(1.95), # treatment effect
                          OR_cond_batchid = 1.25) {
  
  load(rdata_path)
  otu_original <- t(otu)
  p <- ncol(otu_original)

  set.seed(1)
  selected_samples <- sample.int(nrow(otu_original), n)
  
  ## SL 12.2: make treatment and batch correlated via OR
  p_cond  <- 0.5
  p_batch <- 0.5

  # binary correlation
  bin_corr <- bincorr(OR_cond_batchid, p_cond, p_batch)
  cond_batch_mat <- rmvbin(n, margprob = c(p_cond, p_batch), 
    bincorr  = (1 - bin_corr) * diag(2) + bin_corr
  )
  
  #cond <- rbinom(n, 1, 0.5) # treatment/control
  #batchid <- rbinom(n, 1, 0.5) # sex
  cond  <- cond_batch_mat[, 1]
  batchid <- cond_batch_mat[, 2]
  
  ## perturb independent draws of taxa counts for each sample
  otu_tmp <- matrix(nrow = n*m, ncol = p+2)
  otu_tmp <- as.data.frame(otu_tmp)
  names(otu_tmp) <- c('patient.id','time',paste0('taxa',1:p))
  otu_tmp$patient.id <- rep(1:n, each = m)
  otu_tmp$time <- 1:m
  
  ## Initialize taxa counts
  for (i in 1:n) {
    otu_tmp[(i-1)*m + 1, 3:(p+2)] <- bayesm::rdirichlet(otu_original[selected_samples[i],] + 0.5) *
      sum(otu_original[selected_samples[i],])
  }

  ## SL 12.1: compute the prevalence
  otu_mat <- as.matrix(otu_original)
  prevalence <- colSums(otu_mat > 0) / nrow(otu_mat)
  prevalence_order <- order(prevalence, decreasing = TRUE)
  
  ## SL 12.1: assign id_d and id_i based on prevalence order
  id_d_treat <- prevalence_order[1:10]    # first 10 decrease
  id_i_treat <- prevalence_order[11:20]   # last 10 increase
  
  # every other taxa decrease or increase
  # batch_order <- id
  id_d_batch <- prevalence_order[seq(1, length(prevalence_order), by = 2)]
  id_i_batch <- prevalence_order[seq(2, length(prevalence_order), by = 2)]
  
  ## set up the rest of otu table for t=2,3,4
  for (i in 1:n)
  {
    for (j in 2:m)
    {
      # SL 11.14
      prev_row <- (i-1)*m + (j-1)
      cur_row  <- (i-1)*m + j
      
      # perturb previous time points' measurement to create current time point measurements
      prev_counts <- as.numeric(otu_tmp[prev_row, 3:(p+2)])
      otu_tmp[cur_row, 3:(p+2)] <- bayesm::rdirichlet(prev_counts) * sum(prev_counts)

      if (cond[i]) # for treated subjects only
      {
        
        # SL 12.1: update to loop thru each id_i and id_d
        s1_treat <- sum(otu_tmp[cur_row, id_d_treat + 2], na.rm = TRUE)
        s2_treat <- sum(otu_tmp[cur_row, id_i_treat + 2], na.rm = TRUE)
        d <- treat_effect[1]
        
        if (s1_treat > 0 && s2_treat > 0) {
          for (i1 in id_d_treat) {
            otu_tmp[cur_row, i1 + 2] <- otu_tmp[cur_row, i1 + 2] / d
          }
          change_fold <- 1 + (s1_treat/s2_treat)*((d - 1)/d)
          
          for (i2 in id_i_treat) {
            otu_tmp[cur_row, i2 + 2] <- otu_tmp[cur_row, i2 + 2] * change_fold
          }
        }
      }
      
      if (batchid[i]) # for female only
      {
        # SL 12.1: update to loop thru each id_i and id_d
        s1_batch <- sum(otu_tmp[cur_row, id_d_batch + 2], na.rm = TRUE)
        s2_batch <- sum(otu_tmp[cur_row, id_i_batch + 2], na.rm = TRUE)
        d <- batch_effect[1]
        
        if (s1_batch > 0 && s2_batch > 0) {
          for (i1 in id_d_batch) {
            otu_tmp[cur_row, i1 + 2] <- otu_tmp[cur_row, i1 + 2] / d
          }
          change_fold <- 1 + (s1_batch/s2_batch)*((d - 1)/d)
          
          for (i2 in id_i_batch) {
            otu_tmp[cur_row, i2 + 2] <- otu_tmp[cur_row, i2 + 2] * change_fold
          }
        }
      }
    }
  }
  
  otu_tmp_unrounded <- otu_tmp
  # round the otu table
  otu_tmp <- round(otu_tmp)

  example_data <- list(metadata = cbind(otu_tmp[,1:2],
                                        batch = rep(batchid, each = m),
                                        treatment = rep(cond, each = m)) %>%
                         rename(subjectid = `patient.id`) %>%
                         # SL 1.18: 
                         mutate(time = factor(time, 
                                  levels = sort(unique(time)),
                                  labels = paste0("t", sort(unique(time)))
                                              )),
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

simulate_data_sick <- function(rdata_path = "mom_270.Rdata",
                          n = 100, # subjects
                          m = 4, # time point
                          batch_effect = c(1.25), # obscuring sex effect
                          treat_effect = c(1.95), # treatment effect
                          OR_cond_batchid = 1.25) {
  
  load(rdata_path)
  otu_original <- t(otu)
  p <- ncol(otu_original)
  
  set.seed(1)
  selected_samples <- sample.int(nrow(otu_original), n)
  
  ## SL 1.21: Make treatment and batch correlated via OR
  p_cond  <- 0.5
  p_batch <- 0.3
  
  ## Original setup
  # sick_time_options <- expand.grid(rep(list(0:1),m))
  # sick_time <- t(sick_time_options[sample(nrow(sick_time_options), n, replace = TRUE),])
  
  # binary correlation
  bin_corr <- bincorr(OR_cond_batchid, p_cond, p_batch)
  cond_batch_mat <- rmvbin(n, margprob = c(p_cond, p_batch),
                           bincorr  = (1 - bin_corr) * diag(2) + bin_corr
  )
  
  cond  <- cond_batch_mat[, 1]
  sick_base <- cond_batch_mat[, 2]
  
  # make sickness a matrix
  batchid <- matrix(NA, nrow = n, ncol = m)
  
  for (i in 1:n) {
    # baseline sickness from step 1
    base <- sick_base[i]
    
    # make sickness vary across time but keep correlation structure
    batchid[i, ] <- rbinom(m, 1, prob = ifelse(base == 1, 0.7, 0.3))
  }
  batchid <- t(batchid)
  
  ## perturb independent draws of taxa counts for each sample
  otu_tmp <- matrix(nrow = n*m, ncol = p+2)
  otu_tmp <- as.data.frame(otu_tmp)
  names(otu_tmp) <- c('patient.id','time',paste0('taxa',1:p))
  otu_tmp$patient.id <- rep(1:n, each = m)
  otu_tmp$time <- 1:m
  
  ## Initialize taxa counts
  for (i in 1:n) {
    otu_tmp[(i-1)*m + 1, 3:(p+2)] <- bayesm::rdirichlet(otu_original[selected_samples[i],] + 0.5) *
      sum(otu_original[selected_samples[i],])
    
    # sickness initialization
    if (batchid[1,i])
    {
      otu_tmp[(i-1)*m+1, id + 2] <- otu_tmp[(i-1)*m+1, id + 2] * batch_effect[1]
    }
  }
  
  ## SL 12.1: compute the prevalence
  otu_mat <- as.matrix(otu_original)
  prevalence <- colSums(otu_mat > 0) / nrow(otu_mat)
  prevalence_order <- order(prevalence, decreasing = TRUE)
  
  ## SL 12.1: assign id_d and id_i based on prevalence order
  id_d_treat <- prevalence_order[1:10]    # first 10 decrease
  id_i_treat <- prevalence_order[11:20]   # last 10 increase
  
  # every other taxa decrease or increase
  # batch_order <- id
  id_d_batch <- prevalence_order[seq(1, length(prevalence_order), by = 2)]
  id_i_batch <- prevalence_order[seq(2, length(prevalence_order), by = 2)]
  
  ## set up the rest of otu table for t=2,3,4
  for (i in 1:n)
  {
    for (j in 2:m)
    {
      prev_row <- (i-1)*m + (j-1)
      cur_row  <- (i-1)*m + j
      
      # perturb previous time points' measurement to create current time point measurements
      prev_counts <- as.numeric(otu_tmp[prev_row, 3:(p+2)])
      otu_tmp[cur_row, 3:(p+2)] <- bayesm::rdirichlet(prev_counts) * sum(prev_counts)
      
      if (cond[i]) # for treated subjects only
      {
        
        # SL 12.1: update to loop thru each id_i and id_d
        s1_treat <- sum(otu_tmp[cur_row, id_d_treat + 2], na.rm = TRUE)
        s2_treat <- sum(otu_tmp[cur_row, id_i_treat + 2], na.rm = TRUE)
        d <- treat_effect[1]
        
        if (s1_treat > 0 && s2_treat > 0) {
          for (i1 in id_d_treat) {
            otu_tmp[cur_row, i1 + 2] <- otu_tmp[cur_row, i1 + 2] / d
          }
          change_fold <- 1 + (s1_treat/s2_treat)*((d - 1)/d)
          
          for (i2 in id_i_treat) {
            otu_tmp[cur_row, i2 + 2] <- otu_tmp[cur_row, i2 + 2] * change_fold
          }
        }
      }
      
      if (batchid[j, i]) 
      {
        # SL 12.1: update to loop thru each id_i and id_d
        s1_batch <- sum(otu_tmp[cur_row, id_d_batch + 2], na.rm = TRUE)
        s2_batch <- sum(otu_tmp[cur_row, id_i_batch + 2], na.rm = TRUE)
        d <- batch_effect[1]
        
        if (s1_batch > 0 && s2_batch > 0) {
          for (i1 in id_d_batch) {
            otu_tmp[cur_row, i1 + 2] <- otu_tmp[cur_row, i1 + 2] / d
          }
          change_fold <- 1 + (s1_batch/s2_batch)*((d - 1)/d)
          
          for (i2 in id_i_batch) {
            otu_tmp[cur_row, i2 + 2] <- otu_tmp[cur_row, i2 + 2] * change_fold
          }
        }
      }
    }
  }
  
  otu_tmp_unrounded <- otu_tmp
  # round the otu table
  otu_tmp <- round(otu_tmp)
  
  example_data <- list(metadata = cbind(otu_tmp[,1:2],
                                        batch = rep(batchid, each = m),
                                        treatment = rep(cond, each = m)) %>%
                         rename(subjectid = `patient.id`) %>%
                         # SL 1.18: 
                         mutate(time = factor(time, 
                                              levels = sort(unique(time)),
                                              labels = paste0("t", sort(unique(time)))
                         )),
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
