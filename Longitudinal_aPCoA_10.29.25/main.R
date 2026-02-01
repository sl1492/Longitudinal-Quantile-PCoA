source("00_helpers.R")
source("01_simulation.R")
source("02_kernel_pca.R")
source("03_linear_aPCoA.R")
source("04_quantile_aPCoA.R")

# simulate data
sim <- simulate_data_sex(rdata_path = "mom_270.Rdata",
                     n = 100, m = 4,
                     batch_effect = c(1.25), # obscuring sex effect
                     treat_effect = c(1.95), # treatment effect
                     OR_cond_batchid = 1.25)

# build Aitchison kernel & PCs
kpca <- build_kernel_pcs(sim$otu_tmp, sim$example_data)

######### Test #############
PCs_before <- kpca$PCs
var_explained <- apply(PCs_before, 2, var) / sum(apply(PCs_before, 2, var))
round(var_explained[1:5] * 100, 2)
############################

# Linear mixed model
m1 <- lmm_long_apcoa(kpca$PCs, 
                     sim$otu_tmp, 
                     sim$batch, 
                     sim$n, 
                     sim$m)
resid_PCs <- lapply(m1, function(x) process_resids(x, sim$otu_tmp, sim$cond, sim$batch))
p_final1 <- plot_apcoa(resid_PCs[[3]], "LMM")
# Variance explained
resid_PCs[[3]][2:6]
resid_PCs[[3]][[2]]/resid_PCs[[3]][[3]]

# Quantile
m2 <- quantile_apcoa(PCs = kpca$PCs, 
                     metadata = kpca$metadata, 
                     to_remove = c("batch","time"),
                     to_remain = c("treatment"),
                     taus = seq(0.01, 0.99, length.out = 99),
                     tauw = rep(1/99, 99),
                     subject_id = "subjectid"
                     )
final_PCs <- process_resids(m2, sim$otu_tmp, sim$cond, sim$batch)
p_final2 <- plot_apcoa(final_PCs,"Quantile")
# Variance explained
final_PCs[2:6]
final_PCs[[2]]/final_PCs[[3]]

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


