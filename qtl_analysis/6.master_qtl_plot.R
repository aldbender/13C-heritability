library(ggplot2)

data <- read.delim("/5.master_qtl_plot.txt", header=TRUE)

data$trait <- factor(data$trait, levels=c("mindex","pct_iso","pct_anteiso","nalk_ACL","iso_ACL","anteiso_ACL","nalk_CPI","iso_CPI","anteiso_CPI","eps13C_n.C31","eps13C_i.C31","eps13C_n.C33","eps13C_i.C33","n31.i31","n33.i33"))

data$log10_pval <- as.numeric(ifelse(data$p_value>0.05, "NA", -log(data$p_value, base=10)))

data$qual_pval <- ifelse(data$p_value<0.0001, "less_than_0001", ifelse(data$p_value<0.001, "less_than_001", ifelse(data$p_value<0.01, "less_than_01", ifelse(data$p_value<0.05, "less_than_05", "NA"))))

p <- ggplot(data, aes(x=IL, y=trait, fill=qual_pval))
p + geom_tile() + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1)) + scale_fill_manual(values=c("red","orange","yellow","green","gray90"))



