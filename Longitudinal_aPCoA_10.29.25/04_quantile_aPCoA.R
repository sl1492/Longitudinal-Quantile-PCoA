
quantile_apcoa <- function(PCs, metadata,
                           covariates = "all",
                           lambda = 1.5,
                           q = 99,
                           subject_id = "subjectid",
                           formula = "batch + time",  
                           otu_tmp,
                           batchid,
                           cond) {
  
  taus <- seq(0.01, 0.99, length.out = q)
  tauw <- rep(1/q, q)
  df <- cbind(metadata, PCs)
  df[[subject_id]] <- factor(df[[subject_id]])
  
  PC_hat_list <- list()
  
  for (pc in colnames(PCs)) {
    cat("Processing", pc, "...\n")
    
    formula_text <- paste0(pc, " ~ ", formula, " | ", subject_id)
    current_formula <- as.formula(formula_text)
    
    Xmeta <- model.matrix(as.formula(paste0(pc, " ~ ", formula)), data = df)[, -1]
    Xmeta <- as.data.frame(Xmeta)
    y <- df[[pc]]
    
    # SL 10.20: added scale() for covariates
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
    # SL 9.27: Set selected fixed-effect columns to 0 in the design matrix.
    # - If covariates == "all": zero ALL fixed effects (include time) → intercept-only target.
    # - Else: zero only columns whose names match `covariates` and their interactions.
    X_corrected <- X_raw
    
    if (identical(covariates, "all")) {
      
      # use column-wise "null" values (here: minima) for ALL columns
      null_vals <- apply(X_raw, 2, min, na.rm = TRUE)
      X_corrected[,] <- matrix(rep(null_vals, each = nrow(X_raw)), nrow = nrow(X_raw))
    } else {
      
      idx <- grep(paste(covariates, collapse = "|"), colnames(X_raw))
      if (length(idx) > 0) {
        
        # nulls only for selected columns
        null_vals <- apply(X_raw[, idx, drop = FALSE], 2, min, na.rm = TRUE)
        X_corrected[, idx] <- matrix(rep(null_vals, each = nrow(X_raw)), nrow = nrow(X_raw))
      }
    }
    
    Xi_corrected <- cbind(Intercept = 1, as.matrix(X_corrected))
    
    # Get the predicted full quantile curves after correction
    Q_fixed_corrected <- Xi_corrected %*% coef_mat # n × q
    R_subj_corrected  <- random_mat[subject_index, , drop = FALSE]
    quant_matrix_corrected <- Q_fixed_corrected + R_subj_corrected
    colnames(quant_matrix_corrected) <- taus
    
    ### SL 9.29: Quantile matching (invert → transport)
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
  
  final_PCs <- process_resids(resid_PCs, otu_tmp, cond, batchid)
}
