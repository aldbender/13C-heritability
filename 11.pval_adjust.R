library(ggplot2)


#data <- read.table("./10.pvals.txt", header=TRUE)
data <- read.delim("~/GoogleDrive/Grad_school/Research/Tomatoes/Paper_drafts/d13C_CSIA/Figures/3.QTL analysis/QTL modeling/10.pvals.txt")

names(data)

data$bh <- p.adjust(data$pvalue, method="BH")

sig_pvals <- subset(data, data$bh<0.05)

tail(sig_pvals)

write.table(sig_pvals, "12.sig_bh_pvals.txt")


	
	
	