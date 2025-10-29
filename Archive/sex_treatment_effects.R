library(bindata)
library(compositions)
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

D2K <- function(D)
{
  n <- nrow(D)
  centerM <- diag(n) - 1/n
  K <- -0.5 * centerM %*% (D * D) %*% centerM
  eK <- eigen(K, symmetric = TRUE)
  K <- eK$vector %*% diag(pmax(0, eK$values)) %*% t(eK$vector)
  return(K)
}

process_resids <- function(resids, otu_tmp, cond, batchid)
{
  eigen <- eigen(tcrossprod(resids), symmetric = TRUE)
  mK <- which(eigen$values > 1e-9)
  PCs <- eigen$vectors[,mK] %*% diag(sqrt(eigen$values[mK]))/1
  PCs <- as.data.frame(cbind(patient.id = otu_tmp$patient.id,
                             time = otu_tmp$time,
                             treatment = rep(cond, each = max(otu_tmp$time)),
                             batch = rep(batchid, each = max(otu_tmp$time)),
                             PCs))
  names(PCs)[5:ncol(PCs)] <- c(paste0('PC',1:(ncol(PCs)-4)))
  return(list(PCs = PCs,
              PC1_perc = eigen$values[1]/sum(eigen$values[mK]),
              PC2_perc = eigen$values[2]/sum(eigen$values[mK])))
}

load("mom_270.Rdata")
otu_original <- t(otu)
#libsize <- rowSums(otu_original)
p <- ncol(otu_original)
n = 100; m = 4

# change condition taxa
set.seed(1)
cond_taxa1 <- sample(setdiff(1:p, id), 20) # introduce changes in the community profiles shared by all subjects
cond_taxa2 <- sample(setdiff(1:p, cond_taxa1), 20) # Treatment group, could have both sex and treatment effects

set.seed(7711)
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

cond_effect <- c(0.25)
batch_effect <- c(200)
sub_fc <- ifelse(cond == 1, 6*rlnorm(n, meanlog = log(5), sdlog = 1.2), 1) # mean=31.5, median=12.1, max=558.5

for (i in 1:n)
{
  for (j in 2:m)
  {
    # perturb previous time points' measurement to create current time point measurements
    otu_tmp[((i-1)*m + j), 3:(p+2)] <- bayesm::rdirichlet(as.numeric(otu_tmp[((i-1)*m + j - 1), 3:(p+2)] + 0.5)) *
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
save(example_data, file = "example_data.RData")

# Aitchison
rel_tmp <- otu_tmp[,3:ncol(otu_tmp)]/rowSums(otu_tmp[,3:ncol(otu_tmp)])
rel_tmp_clr <- clr_lite(otu_tmp[,3:ncol(otu_tmp)], samples_are = "rows", 
                        method = "logunif", replicates = 100)

K <- D2K(as.matrix(dist(rel_tmp_clr)))
K_eigen <- eigen(K, symmetric = T)
mK <- which(K_eigen$values > 1e-9)
PCs <- K_eigen$vectors[ ,mK] %*% diag(sqrt(K_eigen$values[mK]))


annotations <- data.frame(
  time = c(1, 2, 3, 4),  # Facet labels
  label = c("t[0]", "t[1]", "t[2]", "t[3]"),  # Custom labels
  x = c(-Inf, -Inf, -Inf, -Inf),  # x position (left side)
  y = c(Inf, Inf, Inf, Inf)  # y position (top side)
)

############### linear regression ########################
mK90 <- min(which(cumsum(K_eigen$values[mK]/sum(K_eigen$values[mK])) > 0.9))
PCs <- K_eigen$vectors[ ,1:mK90] %*% diag(sqrt(K_eigen$values[1:mK90]))
stand.marg.resids <- marg.resids <- cond.resids <- matrix(NA, nrow = nrow(PCs), 
                                                          ncol = ncol(PCs))
for (j in 1:ncol(PCs))
{
  dat <- data.frame(y = PCs[,j],
                    time = otu_tmp$time,
                    batch = rep(batchid, each = m),
                    id = rep(1:n, each = m))
  model1 <- lmer(y ~ batch*time + (1 | id), dat)
  if (isSingular(model1)) next
  var.d <- crossprod(getME(model1,"Lambdat")) # relative random-effects covariance-covariance matrix G
  Zt <- getME(model1, "Zt") # random effects design
  vr <- sigma(model1)^2 # estimated residual variance
  var.b <- vr * (crossprod(Zt, var.d) %*% Zt) # Random-effects contribution to the marginal variance ZtGZt^T
  sI <- vr * Diagonal(nrow(dat)) # Noise variance R_i = sigma2I
  var.y <- var.b + sI # Full marginal covariance
  cholVarInv <- solve(t(chol(var.y)))
  
  # Step 5: obtain standardized residuals
  marg.resids[, j] <- dat$y - model.matrix(model1) %*% summary(model1)$coefficients[,1]
  stand.marg.resids[, j] <- as.numeric(cholVarInv %*% marg.resids[, j])
  cond.resids[, j] <- unname(residuals(model1))
}

resid_PCs <- lapply(list(marg.resids, cond.resids, stand.marg.resids),
                    function(x) process_resids(x, otu_tmp, cond, batchid))

p_final1 <- resid_PCs[[3]]$PCs %>%
  mutate(Sex = ifelse(batch == 1, "Female", "Male")) %>%
  mutate(Arm = factor(ifelse(treatment == 0,'Control',"Treatment"), levels = c("Treatment","Control"))) %>%
  ggplot(aes(PC1, PC2, color = Arm, shape = `Sex`, group = `Sex`)) +
  geom_point(size = 2) +
  facet_wrap(~time, nrow = 2, scales = 'free',
             labeller = label_bquote(t[.(time-1)])) +
  theme(strip.text = element_blank()) +
  xlab(paste0('PC1 (', round(100*resid_PCs[[3]]$PC1_perc,2),'%)')) +
  ylab(paste0('PC2 (', round(100*resid_PCs[[3]]$PC2_perc,2),'%)')) +
  theme(plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_blank(), #transparent legend bg
        legend.box.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18)) +
  ggtitle('Standardized Marginal Residual PCA') +
  stat_ellipse(aes(group = `Sex`, lty = `Sex`)) +
  geom_text(data = annotations, aes(x = x, y = y, label = label),
            hjust = -0.2, vjust = 1.5, parse = TRUE, inherit.aes = FALSE)
p_final1


############### quantile regression ########################
colnames(PCs)[1:mK90] <- paste0("PC", 1:mK90)
df <- cbind(example_data$metadata, PCs)

q <- 299
taus <- ppoints(q + 2)[2:(q + 1)]
tauw <- rep(1/q, q)            # or a mild upper-tail weight
lambda <- 0.25
#q <- 99 # quantiles
#lambda <- 10 # penalty parameter controlling the shrinkage
#tauw <- rep(1/q, q)
#taus <- seq(0.01, 0.99, length.out = q)

formula <- "batch + time + batch*time"
subject_id <- "subjectid"
covariates <- c("all")
df$subjectid <- factor(df$subjectid)
quantiles <- seq(0.01, 0.99, length.out = q)

PC_hat_list <- list()

for (pc in colnames(PCs)) {
  cat("Processing", pc, "...\n")

  formula_text <- paste0(pc, " ~ ", formula, " | ", subject_id)
  current_formula <- as.formula(formula_text)
  
  metadata <- model.matrix(as.formula(paste0(pc, " ~ ", formula)), data = df)[,-1]
  metadata <- as.data.frame(metadata)
  y <- df[[pc]]
  X_raw <- model.matrix(~., metadata)[,-1]
  
  # scaling deign matrix
  #X.sc <- scale(X_raw, center=T, scale=T)
  #X.sc.q = X.sc[, colSums(!is.na(X.sc)) > 0, drop=F]
  X.sc.q = X_raw
  
  len_cov <- ncol(X.sc.q) + 1
  
  # Fit rqpd
  fit <- rqpd(formula = current_formula, data = df,
              control = quantreg::sfn.control(tmpmax = 1e8),
              panel(lambda = lambda,
                    tauw = tauw,
                    taus = taus)
  )
  #coefs <- fit$coefficients
  # Extract coefficient
  coef_fixed = fit$coef[1:(len_cov*q)]
  coef_mat = matrix(coef_fixed, nrow=len_cov, ncol=q, byrow=F)
  
  coef_random <- fit$coef[-(1:(len_cov * q))] 
  subject_levels <- levels(df[[subject_id]])  # or however subjects were coded
  subject_index <- as.numeric(factor(df[[subject_id]], levels = subject_levels))
  random_mat <- matrix(coef_random, nrow = length(subject_levels), ncol = q, byrow = FALSE)
  
  # Predict quantiles (with batch)
  quant_matrix <- matrix(NA, nrow = nrow(X.sc.q), ncol = q)
  for (i in 1:nrow(X.sc.q)) {
    xi <- as.numeric(c(1, X.sc.q[i, ]))
    subj_idx <- subject_index[i]  # gives 1 to 100
    ri <- random_mat[subj_idx, ]  # 1 × q
    quant_matrix[i, ] <- xi %*% coef_mat + ri
  }
  
  # Predict quantiles (batch free)
  X_corrected <- X.sc.q
  if (identical(covariates, "all")) {
    X_corrected[,] <- 0
  } else {
    batch_covariate_idx <- grep(paste(covariates, collapse = "|"), colnames(X.sc.q))
    X_corrected[, batch_covariate_idx] <- 0
  }
  colnames(quant_matrix) <- quantiles
  
  quant_matrix_corrected <- matrix(NA, nrow = nrow(X.sc.q), ncol = q)
  for (i in 1:nrow(X_corrected)) {
    xi <- c(1, X_corrected[i, ])
    subj_idx <- subject_index[i]  # gives 1 to 100
    ri <- random_mat[subj_idx, ]  # 1 × q
    quant_matrix_corrected[i, ] <- xi %*% coef_mat + ri
  }
  colnames(quant_matrix_corrected) <- quantiles
  
  PC_hat <- vapply(seq_along(y), function(i) {
    q_obs <- quant_matrix[i, ]
    q_corr <- quant_matrix_corrected[i, ]
    
    # Find quantiles where predicted value is ≤ observed
    match_idx <- which(q_obs <= y[i])
    
    if (length(match_idx) > 0) {
      tau_idx <- max(match_idx)  # left-continuous quantile
    } else {
      tau_idx <- 1  # fallback to smallest τ
    }
    q_corr[tau_idx]
  }, numeric(1))
  PC_hat_list[[pc]] <- PC_hat
}

resid_PCs <- do.call(cbind, PC_hat_list)
colnames(resid_PCs) <- colnames(PCs)
final_PCs <- process_resids(resid_PCs, otu_tmp, cond, batchid)


############# Plot 2 ############
p_final2 <- final_PCs$PCs %>%
  mutate(Sex = ifelse(batch == 1, "Female", "Male")) %>%
  mutate(Arm = factor(ifelse(treatment == 0,'Control',"Treatment"), levels = c("Treatment","Control"))) %>%
  ggplot(aes(PC1, PC2, color = Arm, shape = `Sex`, group = `Sex`)) +
  geom_point(size = 2) +
  facet_wrap(~time, nrow = 2, scales = 'free',
             labeller = label_bquote(t[.(time-1)])) +
  theme(strip.text = element_blank()) +
  xlab(paste0('PC1 (', round(100*final_PCs$PC1_perc,2),'%)')) +
  ylab(paste0('PC2 (', round(100*final_PCs$PC2_perc,2),'%)')) +
  theme(plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_blank(), #transparent legend bg
        legend.box.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18)) +
  ggtitle('Quantile Adjusted Longitudinal PCoA ') +
  stat_ellipse(aes(group = `Sex`, lty = `Sex`)) +
  geom_text(data = annotations, aes(x = x, y = y, label = label),
            hjust = -0.2, vjust = 1.5, parse = TRUE, inherit.aes = FALSE)
p_final2



############ final plot #################
# Plot for extracting the 'Sex' legend
p_final1_sick <- p_final1 +
  scale_color_manual(values = c("Treatment" = "#F8766D", "Control" = "#00BFC4")) +  # Correct color values
  guides(color = "none") +  # Hide the color legend
  guides(shape = guide_legend(title = "Sick"))  # Only show shape legend

# Plot for extracting the 'Arm' legend
p_final1_arm <- p_final1 +
  scale_shape_manual(values = c("True" = 16, "False" = 17)) +  # Correct shape values
  guides(shape = "none", lty = "none") +  # Hide the shape and line type (lty) legend
  guides(color = guide_legend(title = "Arm"))  # Only show color legend

legend_sick <- get_legend(p_final1_sick)
legend_arm <- get_legend(p_final1_arm)
combined_legend <- plot_grid(legend_sick, legend_arm, ncol = 2, align = "v")
p_final1_edited <- p_final1 +
  theme(legend.position = "none",
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  ggtitle("")
p_final2_edited <- p_final2 +
  theme(legend.position = "none",
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  ggtitle("")

final_plot <- plot_grid(
  plot_grid(p_final2_edited, p_final1_edited, ncol = 1, align = 'v', 
            labels = c("(a)", "(b)"), label_size = 18),
  combined_legend,
  ncol = 1,
  rel_heights = c(1, 0.15)  # Adjust height ratio between plot and legend
)
ggdraw(final_plot)

legend <- get_legend(p_final1)
# Stack quantile and linear plots vertically
plots_stacked <- plot_grid(
  p_final2_edited,
  p_final1_edited,
  ncol = 1,
  labels = c("quantile", "linear")
)

# Combine stacked plots and legend side by side
ggdraw(plot_grid(
  plots_stacked,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.25)  # Adjust width of plot vs legend area
))

