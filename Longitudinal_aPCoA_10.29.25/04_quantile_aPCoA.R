quantile_apcoa <- function(PCs, metadata,
                           to_remove = c("batch","time"),
                           to_remain = c("treatment"),
                           q = 99,
                           taus = seq(0.01, 0.99, length.out = q),
                           tauw = rep(1/q, q),
                           subject_id = "subjectid"
                           ) {
  
  df <- cbind(metadata, PCs)
  df[[subject_id]] <- factor(df[[subject_id]])
  
  PC_hat_list <- list()
  
  for (pc in colnames(PCs)) {
    cat("Processing", pc, "...\n")
    
    ## SL 1.5: Compute lambda for each PC
    dat <- df[pc]
    dat_scale <- scale(dat, center = T, scale = T)
    # Group by subjectid
    dat_mean <- aggregate(dat_scale, by = list(subjectid = df[,"subjectid"]), FUN = mean)[,-1]
    lambda <- 1 / sd(dat_mean, na.rm = TRUE)
    
    ## SL 1.18: Reconstruct formula using re_remain and to_remove
    formula <- paste(c(to_remove, to_remain), collapse = " + ")
    formula_text <- paste0(pc, " ~ ", formula, " | ", subject_id)
    current_formula <- as.formula(formula_text)
    
    Xmeta <- model.matrix(as.formula(paste0(pc, " ~ ", formula)), data = df)[, -1]
    Xmeta <- as.data.frame(Xmeta)
    y <- df[[pc]]
    
    X_raw <- scale(model.matrix(~., Xmeta)[,-1])
    df[, colnames(X_raw)] <- X_raw
    
    len_cov <- ncol(X_raw) + 1
    
    ### Fit rqpd
    fit <- rqpd(formula = current_formula, data = df,
                control = quantreg::sfn.control(tmpmax = 1e8),
                panel(lambda = lambda,
                      tauw = tauw,
                      taus = taus))
    
    # Extract coefficient
    coef_fixed = fit$coef[1:(len_cov*q)]
    coef_mat = matrix(coef_fixed, nrow=len_cov, ncol=q, byrow=F)
    
    coef_random <- fit$coef[-(1:(len_cov * q))]
    subject_levels <- levels(df[[subject_id]])  
    subject_index <- as.numeric(factor(df[[subject_id]], levels = subject_levels))
    random_mat <- matrix(coef_random, nrow = length(subject_levels), ncol = q, byrow = FALSE)
    
    ### Predict quantiles
    ## Predict quantiles (with batch)
    
    Xi <- cbind(Intercept = 1, as.matrix(X_raw))
    # Fixed-effects
    Q_fixed <- Xi %*% coef_mat
    # Subject random effects
    R_subj  <- random_mat[subject_index, , drop = FALSE]
    # Full fitted conditional quantile curves
    quant_matrix <- Q_fixed + R_subj
    colnames(quant_matrix) <- taus
    
    ## Predict quantiles (Remove target covariate)
    X_corrected <- X_raw
    
    # SL 1.18: Replace selected covariates with their minima, keeping treatment (to_remain)
    if (length(to_remove) > 0) {
      # find column indices matching the to_remove covariates
      idx_remove <- grep(paste(to_remove, collapse = "|"), colnames(X_raw), ignore.case = TRUE)
      
      if (length(idx_remove) > 0) {
        # replace those columns in X_corrected with their minima        
        null_vals <- apply(X_raw[, idx_remove, drop = FALSE], 2, min, na.rm = TRUE)
        X_corrected[, idx_remove] <- matrix(rep(null_vals, each = nrow(X_raw)),
                                            nrow = nrow(X_raw))
      }
    }
    
    # Keep to_remain covariates unchanged
    if (length(to_remain) > 0) {
      idx_remain <- grep(paste(to_remain, collapse = "|"), colnames(X_raw), ignore.case = TRUE)
      if (length(idx_remain) > 0) {
        X_corrected[, idx_remain] <- X_raw[, idx_remain, drop = FALSE]
      }
    }
    
    # add intercept for fixed effects
    Xi_corrected <- cbind(Intercept = 1, as.matrix(X_corrected))
    
    # Get the predicted full quantile curves after correction
    Q_fixed_corrected <- Xi_corrected %*% coef_mat # n × q
    R_subj_corrected  <- random_mat[subject_index, , drop = FALSE]
    quant_matrix_corrected <- Q_fixed_corrected + R_subj_corrected
    colnames(quant_matrix_corrected) <- taus
    
    ### Quantile matching (invert → transport)
    PC_hat <- vapply(seq_along(y), function(i) {
      
      # Fitted conditional quantile curve for observation i
      q_obs <- quant_matrix[i, ]         
      
      # Corrected quantile curve for the same observation
      q_corr <- quant_matrix_corrected[i, ]
      
      ## Invert step
      # Find all τ grid points where the fitted quantile is ≤ the observed outcome y[i].
      match_idx <- which(q_obs <= y[i])
      
      ## Left-continuous inverse:
      # take the largest τ index such that Q̂(τ | x_i) ≤ y[i].
      # If none exist (y[i] smaller than all fitted values), fall back to τ_1.
      if (length(match_idx) > 0) {
        tau_idx <- max(match_idx) 
      } else {
        tau_idx <- 1
      }
      
      ## Transport step
      # Map the same conditional rank τ̂_i to the target quantile curve.
      # This gives the batch-corrected value for observation i.
      q_corr[tau_idx]
      
    }, numeric(1)) 
    
    PC_hat_list[[pc]] <- PC_hat
  }
  
  resid_PCs <- do.call(cbind, PC_hat_list)
  colnames(resid_PCs) <- colnames(PCs)
  
  resid_PCs
}
