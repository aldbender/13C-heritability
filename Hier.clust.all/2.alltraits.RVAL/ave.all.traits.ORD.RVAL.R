#Load ggplot2
library(ggplot2)

#IMPORTANT: The goal of this script is to generate a heatmap of the rho values from the previous script in which the traits are ordered in the same order as the hierarchical clustering dendrogram. To do that, you will have to copy and paste the traits from the dendrogram into excel and assign them an order number. Once they have an order number, each trait pair, Trait1 and Trait2, for which there is a rho value, must have the associated order attached to it. Please look at the file that is read in below to get a feel for what I mean by this. You will have to use the merge() function to merge the orders for each trait pair correctly

#Read in data
data <- read.table("/traits.ORD.RVAL.txt",header=TRUE)
names(data)

#Make a heat map
t <- ggplot(data=data, aes(x=ord.trait1, y=ord.trait2, fill=r.value))

t + geom_tile() + scale_fill_gradient2(low="firebrick", mid="white", high="forestgreen")  +theme_bw() + theme(axis.text.x=element_text(angle=-90)) + coord_fixed()
t + geom_tile() + scale_fill_gradient2(low="red", mid="black", high="yellow")  +theme_bw() + theme(axis.text.x=element_text(angle=-90)) + coord_fixed()

#NOTE: If you would like to use different colors, check out colorbrewer.org and copy and paste in the hex code colors

# Below: plot without anything on the edges and no legend!
t + geom_tile() + scale_fill_gradient2(low="firebrick", mid="white", high="forestgreen")  +theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + coord_fixed() + guides(fill=FALSE)



