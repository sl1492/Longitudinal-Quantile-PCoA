library(phyloseq)
library(HMP2Data)

mom = momspi16S()
mom = prune_samples(momspi16S_samp$sample_body_site=="vagina" &
                    momspi16S_samp$visit_number==4,
                    mom)
# prune to samples with at least 4000 bugs
mom = prune_samples(sample_sums(mom)>=4000, mom)
# taxa zero rate
zp = rowMeans(otu_table(mom)==0)
# only keep taxa such that at least zero rate is less than 95%
mom = prune_taxa(zp <= 0.95, mom)

otu = matrix(otu_table(mom), ncol=nsamples(mom))
rownames(otu) = taxa_names(mom)
colnames(otu) = sample_names(mom)





id = NULL
zp = rowMeans(otu==0)

# get the 2 taxa that have zero rate...
# ...closest to 0, 0.1, 0.2,...
for (tau in 0:9/10) {
  id = c(id,order(abs(zp-tau))[1:2])
}

temp = zp[id]
names(temp) = paste("taxa", id)
names(temp) = NULL
pdf("/Users/amarise/Documents/Longitudinal Microbiome Visualization/sims from Jiuyao/zerorate_differential.pdf", width=7.5, height=5)
par(font.main=1)
barplot(temp, ylim=c(0, 1), main="Zero Rate of 20 Chosen Differentially Abundant Taxa")
dev.off()

# estimate probabilities from DM distr
shape = dirmult::dirmult(t(otu))$gamma
# num of bugs for each sample
libsize = as.vector(colSums(otu))


id1 = 1:5 # 5 taxa - low zero rates (0-0.2)
id2 = 6:20 # 15 taxa - zero rates 0.2-0.9
id3 = 11:15 # 5 taxa - 0.5-0.7
id4 = 16:20 # 5 taxa - 0.7-0.9

# id consists of reprentative taxa across differential zero rates

save(mom, otu, id,
     shape, libsize,
     id1,id2,id3,id4,
     file="/Users/amarise/Documents/Longitudinal Microbiome Visualization/sims from Jiuyao/mom_270.Rdata")
