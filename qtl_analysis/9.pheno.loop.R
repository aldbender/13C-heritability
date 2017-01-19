library(ggplot2)


#data <- read.table("./8.avg.rpm.values.txt", header=TRUE)
data <- read.delim("~/GoogleDrive/Grad_school/Research/Tomatoes/Paper_drafts/d13C_CSIA/Figures/3.QTL analysis/QTL modeling/8.avg.rpm.values.txt")

names(data)
attach(data)
summary(data)
head(data)

tdata <- as.data.frame(t(data[2:77]))
colnames(tdata) <- t(data[,1])
names(tdata)


###### FOR  REAL#########

#overall.table <- matrix(nrow=8*20332, ncol=4)
overall.table <- matrix(nrow=16*20332, ncol=4)

for(j in c(1:20332 ) ) {
	
#quant.table <- matrix(nrow=8, ncol=4)
  quant.table <- matrix(nrow=16, ncol=4)
  
#for (i in c(20333:20340) ) {
for (i in c(20333:20348) ) {

	print(j)

	t <- lm((tdata[,j]) ~ tdata[,i])
	slope <- summary(t)[[4]][[2]]
	pvalue <- summary(t)[[4]][[8]]

	quant.table[i-20332,1] <- names(tdata)[j]
	quant.table[i-20332,2] <- names(tdata)[i]
	quant.table[i-20332,3] <- slope
	quant.table[i-20332,4] <- pvalue
	
	}
	
	#begin.row <- 1 + (j-1)*8
	#end.row <- begin.row + 7
  begin.row <- 1 + (j-1)*16
  end.row <- begin.row + 15

	overall.table[begin.row:end.row,] <- quant.table
	
	}
	

	
colnames(overall.table)	<- c("itag","cor.trait", "slope","pvalue")

head(overall.table)
tail(overall.table)

write.table(overall.table, "10.pvals.txt")




	
	
	
