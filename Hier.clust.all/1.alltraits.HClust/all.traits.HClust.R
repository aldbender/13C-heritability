library(stats)
library(geneplotter)
library(gplots)
library(Hmisc)
library(reshape)
library(ggplot2)

#Read in data, look at column names with head() and names()
data <- read.table("~/1.ave.alltraits.txt",header=TRUE)
head(data)
names(data)

#Select just those columns of traits that you would like to create a hierarchical clustering dengrogram for
n.data <- data[,c(2:238)]

#Rename the rows with the il number
row.names(n.data) <- data$il

#Convert the dataframe into a matrix to calculate p-values
mn.data <- as.matrix(n.data)

#Use the rcorr function to calculate an r (or in this case Spearman's rho, which is always safer but more conservative)
a <- rcorr(mn.data, type="spearman")

#You will see that you have created a matrix table of rho, p values, and other info
a

#But we need this in a column form, with one column for one trait, the other column the other trait, and the p value for the correaltion in the other column. Using the melt function on $P (the p value) we can achieve this
melt.a <- melt(a$P)

#Now, we have p values, but we need them multiple-test adjusted. Use the p.adjust function to calculate "BH" adjusted p values
BH <- p.adjust(melt.a$value, method=c("BH"))

#Great, now let's write out and save our p values and BH-adjusted p values
write.table(cbind(melt.a,BH), file="all.traits.PVAL.txt")

#Just as we have written out and saved our p values, we need to do the same for our r values (in this case, rho)
melt.b <- melt(a$r)
write.table(melt.b, file="all.traits.RVAL.txt")

#Finally, create and save a hierarchical cluatering dendrogram
snd <- scale(n.data)

hc <- hclust(as.dist(1-abs(cor(snd, method="spearman", use="pairwise.complete.obs"))), method="ward")

plot(hc, cex=1)

#IMPORTANT: Do save the dendrogram for the figure itself. BUT THE MOST IMPORTANT THING WE NEED IS THE ORDER OF THE TRAITS, LEFT TO RIGHT, ALONG THE BOTTOM OF THE CLUSTERING. Save the PDF as large as possible (that is, blow it up to a large size on your screen before saving it). Unfortunately I have never found a good way using code to get the order of the traits. What I do is save a PDF, open up the PDF in Acrobat, and then using the text cursor, highlight the text and copy it from the PDF and then paste in Excel and get the order from the rows. Once you have that order, each trait will have a numbered ordered. For Trait1 and Trait2 and associated p values and r/rho values, then each will have an order which you can use to make the heatmap.
