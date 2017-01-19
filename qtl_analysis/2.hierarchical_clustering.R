library(stats)
library(geneplotter)
library(gplots)
library(Hmisc)
library(reshape)
library(ggplot2)
library(ape)
library(ggdendro)

data <- read.delim("/01.master_qtl_waxsummary.txt",header=TRUE)
head(data)
names(data)


####Generating hierarchical clustering####

n.data <- data[,c(2:17)]

snd <- scale(n.data)

hc <- hclust(as.dist(1-abs(cor(snd, method="spearman", use="pairwise.complete.obs"))), method="ward")

plot(hc, cex=0.8)

plot(as.phylo(hc), type="cladogram", cex=1.2, label.offset=0.01)


###GENERATRE PVALUE MATRIX###
mn.data <- as.matrix(n.data)

a <- rcorr(mn.data, type="spearman")
melt.a <- melt(a$P)
BH <- p.adjust(melt.a$value, method=c("BH"))
pval <- cbind(melt.a,BH)

rho <- melt(a$r)

c <- cbind(pval, rho)

write.table(c, "4.rho_pvals.txt")

