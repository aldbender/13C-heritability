#Load ggplot2
library(ggplot2)

#IMPORTANT: The goal of this script is to generate a heatmap of the p values values from the previous script in which the traits are ordered in the same order as the hierarchical clustering dendrogram. To do that, you will have to copy and paste the traits from the dendrogram into excel and assign them an order number. Once they have an order number, each trait pair, Trait1 and Trait2, for which there is a rho value, must have the associated order attached to it. Please look at the file that is read in below to get a feel for what I mean by this. You will have to use the merge() function to merge the orders for each trait pair correctly

#Also note that for this script, unlike the one for rho, that there are additional columns, "BH" and "BH.NA". The "NA"s are added replaced for those BH p values > 0.05, the siginificance value. Also note that BH p values that are 0 are convert to 2.2x10^-16, the minimum value that R calculates. This will get around the problem of log10 transforming the p value below

#Read in data
data <- read.table("~/traits.ORD.PVAL.txt",header=TRUE)
names(data)

#Make a heat map. Note that the fill is -log10 of BH.NA. The "NA"s will not be color coded, because they are not significant
  
t <- ggplot(data=data, aes(x=ord.trait1, y=ord.trait2, fill=-log(BH.NA,10)))

t + geom_tile() + scale_fill_gradient2(low="darkorange4", mid="orange", high="purple4")  +theme_bw() + theme(axis.text.x=element_text(angle=-90)) + coord_fixed()

#NOTE: If you would like to  use different colors, check out colorbrewer.org and copy and paste in the hex code colors


