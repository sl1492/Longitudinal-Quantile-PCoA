lmm_long_apcoa <- function(PCs, 
                           otu_tmp, 
                           batch, 
                           n, m) {
  
  stand.marg.resids <- marg.resids <- cond.resids <- matrix(NA, nrow = nrow(PCs), ncol = ncol(PCs))

  for (j in 1:ncol(PCs)) {
    
    dat <- data.frame(y = PCs[,j],
                      time = otu_tmp$time,
                      batch = batch,
                      id = rep(1:n, each = m))
    ### fit the linear model
    # using dummy time
    model1 <- lmer(y ~ batch + time + (1 | id), dat)
    
    if (isSingular(model1)) next
    var.d <- crossprod(getME(model1,"Lambdat")) # relative random-effects covariance-covariance matrix G
    Zt <- getME(model1, "Zt") # random effects design
    vr <- sigma(model1)^2 # estimated residual variance
    var.b <- vr * (crossprod(Zt, var.d) %*% Zt) # Random-effects contribution to the marginal variance ZtGZt^T
    sI <- vr * Diagonal(nrow(dat)) # Noise variance R_i = sigma2I
    var.y <- var.b + sI # Full marginal covariance
    cholVarInv <- solve(t(chol(var.y)))
    
    ### Obtain standardized residuals
    marg.resids[, j] <- dat$y - model.matrix(model1) %*% summary(model1)$coefficients[,1]
    stand.marg.resids[, j] <- as.numeric(cholVarInv %*% marg.resids[, j])
    cond.resids[, j] <- unname(residuals(model1))
  }
  
  list(
    marginal = marg.resids,
    conditional = cond.resids,
    standardized = stand.marg.resids
  )
}
